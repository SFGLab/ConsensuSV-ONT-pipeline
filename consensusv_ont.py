import argparse
import vcf
import pysam
import cv2
import numpy as np
import os
import time
import pandas as pd
from tqdm import tqdm
from collections import defaultdict

import tensorflow as tf
from tensorflow.keras import layers, models
from tensorflow.keras.utils import Sequence


def arg_parse():
    """Parsing the arguments.
    Returns:
        argparse.ArgumentParser: Argument parser from argparse.
    """
    parser = argparse.ArgumentParser(description='Gets the SV consensus.')

    parser.add_argument('-i', '--input', help='Input .vcf file.', required=True)
    parser.add_argument('-b', '--bam', help='Input .bam file.', required=True)
    parser.add_argument('-o', '--output', help='Output file name.', default="consensuSV-ONT.vcf")
    parser.add_argument('-t', '--type', help='Mutation type. DEL or INS', default="DEL")        
    parser.add_argument('-s', '--sample', help='Sample name', required=True)
    
    args = parser.parse_args()
    return args


class DataGenerator(Sequence):
    def __init__(self, x_set, y_set, batch_size):
        self.x, self.y = x_set, y_set
        self.batch_size = batch_size

    def __len__(self):
        return int(np.ceil(len(self.x) / float(self.batch_size)))

    def __getitem__(self, idx):
        batch_x = self.x[idx * self.batch_size:(idx + 1) * self.batch_size]
        batch_y = self.y[idx * self.batch_size:(idx + 1) * self.batch_size]
        return batch_x, batch_y

class ImageEncoder:
    
    def encode_to_image(self, row_id, image, image_start, start, end, sv_len, encoded_row):
        if image_start >= start: # if read starts before start of 3 x sv_len window
            encode_start = image_start-start
            encode_end = image_start-start+sv_len*3
            encoded_image = encoded_row[encode_start:encode_end]
            image[row_id, 0:encoded_image.shape[0], :] = encoded_image
        
        else:
            encode_start = abs(image_start-start)
            encode_end = sv_len*3 - encode_start
            encoded_image = encoded_row[0:encode_end]
            image[row_id, encode_start:encode_start+encoded_image.shape[0], :] = encoded_image
        
        return image
    
    def encode(self, chr_, bp, bam_file, mutation_type="DEL"):
        sv_len = bp[1] - bp[0]

        rows_dict = defaultdict(list)
        
        #iterate over reads mapped to mutation region
        for read in bam_file.fetch(chr_, bp[0], bp[1], multiple_iterators=True):
            # filterout low quality reads
            if read.mapping_quality < 60:
                continue
            
            counter = 0
            
            #initialize empty array (zeros - black color) depends of read lenght
            row = np.zeros((read.reference_end-read.reference_start, 3))
            
            for ct in read.cigartuples:
                # Red color as a mapped nucleotides
                if ct[0] in [0]:
                    row[counter:counter+ct[1], :] = [255, 0, 0]
                    counter += ct[1]
                # Green color as a deletions
                elif ct[0] in [2]:
                    row[counter:counter+ct[1], :] = [0, 255, 0]
                    counter += ct[1]
                # Blue color as a insterions. We encode insertions only for mutations_type=INS and only in mutation region, not in extended space
                elif ct[0] in [1] and (counter >= bp[0] - read.reference_start - sv_len and counter <= bp[1] - read.reference_start + sv_len) and mutation_type == "INS":
                    row[counter:counter+ct[1], :] = [0, 0, 255]
                    counter += ct[1]

            rows_dict[read.query_name].append([row, read.reference_start,read.reference_end])

        if not rows_dict:
            return None
        #print(len(rows_dict), bp)
        #initialize variant image with zeros
        image = np.zeros(
            (
                len(rows_dict.keys()), 
                sv_len*3, 
                3
            )
        )
        
        # image start relative to position on chromosome  
        image_start = bp[0] - sv_len

        row_id = 0
        for key, value in rows_dict.items():
            # encode only reads mapped to single region 
            if len(value) == 1:
                encoded_row, start, end = value[0]
                image = self.encode_to_image(row_id, image, image_start, start, end, sv_len, encoded_row)
                row_id += 1
            elif len(value) == 2:
                value = sorted(value, key=lambda x: x[1])
                el1 = value[0]
                el2 = value[1]
                image = self.encode_to_image(row_id, image, image_start, el1[1], el1[2], sv_len, el1[0])
                image = self.encode_to_image(row_id, image, image_start, el2[1], el2[2], sv_len, el2[0])
                image[row_id, el1[2]-image_start:el2[1]-image_start, :] = [0, 255, 0]
                row_id += 1
        # if less than two valid reads mapped - return None
        if row_id < 2:
            return None
        
        resized_image = cv2.resize(image[:row_id, :], (50, 50))
        return np.array(resized_image) /255
    
    
class ConsensuSVONTCore:
    def __init__(self):
        self.encoder = ImageEncoder()
        self.deletion_model_weight = "/tools/ConsensusSV-ONT-pipeline/weights/deletion_model.hdf5"
        self.insertion_model_weight = "/tools/ConsensusSV-ONT-pipeline/weights/insertion_model.hdf5"
        self.valid_chromosomes = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY".split(",")
    
    def get_del_model(self):
        model = models.Sequential()
        model.add(layers.Conv3D(4, 3, activation='relu', padding='same', input_shape=(50, 50, 3, 1)))
        model.add(layers.BatchNormalization())
        model.add(layers.MaxPooling3D((2, 2, 1)))
        model.add(layers.Conv3D(8, (3, 3, 3), activation='relu', padding='same'))
        model.add(layers.BatchNormalization())
        model.add(layers.MaxPooling3D((2, 2, 1)))
        model.add(layers.Conv3D(16, (3, 3, 3), activation='relu', padding='same'))
        model.add(layers.BatchNormalization())
        model.add(layers.MaxPooling3D((2, 2, 1)))
        model.add(layers.Flatten())
        model.add(layers.Dense(32, activation='relu'))
        model.add(layers.Dense(1, activation="sigmoid"))
        optimizer=tf.keras.optimizers.legacy.Adam(learning_rate=0.001)

        model.compile(optimizer=optimizer,
                      loss=tf.keras.losses.BinaryCrossentropy(),
                      metrics=["AUC"],
                     )
        model.load_weights(self.deletion_model_weight)
        return model

    def get_ins_model(self):
        model = models.Sequential()
        model.add(layers.Conv3D(4, 3, activation='relu', padding='same', input_shape=(50, 50, 3, 1)))
        model.add(layers.BatchNormalization())
        model.add(layers.MaxPooling3D((2, 2, 1)))
        model.add(layers.Conv3D(4, (3, 3, 3), activation='relu', padding='same'))
        model.add(layers.BatchNormalization())
        model.add(layers.MaxPooling3D((2, 2, 1)))
        model.add(layers.Flatten())
        model.add(layers.Dense(16, activation='relu'))
        model.add(layers.Dropout(0.2))
        model.add(layers.Dense(1, activation="sigmoid"))

        optimizer=tf.keras.optimizers.legacy.Adam(learning_rate=0.001)
        model.compile(optimizer=optimizer,
                      loss=tf.keras.losses.BinaryCrossentropy(),
                      metrics=["AUC"],
                     )
        model.load_weights(self.insertion_model_weight)
        return model
    
    def read_vcf(self, vcf_file, type_="DEL"):
        vcf_reader = vcf.Reader(filename=vcf_file)
        variants = []
        for record in vcf_reader:
            if record.INFO.get("SVTYPE") == type_ and record.CHROM in self.valid_chromosomes:
                variants.append(record)
        return variants
    
    def write_vcf(self, dataframe, output_filename, sample_id):
        dataframe = dataframe.sort_values(["CHROM", "POS"])
        vcf_header = open("/tools/ConsensusSV-ONT-pipeline/header.txt", "r").read().replace("SAMPLENAME", sample_id)
        
        vcf_lines = []
        for index, row in dataframe.iterrows():
            vcf_lines.append(f"{row['CHROM']}\t{row['POS']}\t.\tN\t<{row['SV_TYPE']}>\t.\tPASS\tSVLEN={row['SVLEN']};SVTYPE={row['SV_TYPE']};END={row['END']};SCORE={row['score']}\tGT\t./.")

        with open(output_filename, 'w') as f:
            f.write(vcf_header)
            f.write('\n'.join(vcf_lines))
    
    def get_model(self, type_):
        if type_ == "DEL":
            return self.get_del_model()
        elif type_ == "INS":
            return self.get_ins_model()
    
    def encode_variants(self, variants, bam_file, mutation_type):
        images, labels = [], []
        for sv in tqdm(variants):
            if isinstance(sv.INFO.get("SVLEN"), list):
                sv_len = abs(sv.INFO.get("SVLEN")[0])
            else:
                sv_len = abs(sv.INFO.get("SVLEN"))
            bp = [sv.POS, sv.POS+sv_len]
            
            # skip huge variants
            if sv_len > 1_000_000 or sv_len < 50:
                continue
            
            img = self.encoder.encode(sv.CHROM, bp, bam_file, mutation_type)
            if img is not None:
                images.append(img)
                labels.append([sv.CHROM] + bp + [sv_len])
        
        return DataGenerator(np.array(images), np.array([0]*len(images)), 32), labels

    def filter_variants(self, vcf_file, bam_path, type_, save_name, sample_id):
        model = self.get_model(type_)
        bam_file = pysam.AlignmentFile(bam_path, "rb")
        variants = self.read_vcf(vcf_file, type_)[:100]
        encoded_variants, encoded_positions = self.encode_variants(variants, bam_file, type_)

        pred_score = model.predict(encoded_variants)

        df_variants = pd.DataFrame(encoded_positions, columns=["CHROM", "POS", "END", "SVLEN"])
        df_variants["SV_TYPE"] = type_
        df_variants["score"] = pred_score
        
        df_variants_filtered = df_variants.loc[df_variants.score >= 0.5]
        
        self.write_vcf(df_variants_filtered, save_name, sample_id)        
        
if __name__ == "__main__":
    args = arg_parse()
    
    s = time.time()
    consensusSV =  ConsensuSVONTCore()
    consensusSV.filter_variants(args.input, args.bam, args.type, args.output, args.sample)
    e = time.time()
    
    print("Processing time:", round((e-s)/60, 2), "minutes") 
