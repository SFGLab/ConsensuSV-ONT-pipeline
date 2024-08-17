# ConsensusSV-ONT-pipeline

## What is ConsensusSV-ONT?
The tool designed for getting consensus out of multiple SV callers' results. The method uses six independent, state-of-the-art structural variant callers for long-read sequencing along with a convolutional neural network for filtering high-quality variants. We provide a runtime environment in the form of a docker image, wrapping a nextflow pipeline for efficient processing using parallel computing.

Docker image: https://hub.docker.com/repository/docker/antonipietryga/consensussv-ont-pipeline/
```
docker pull antonipietryga/consensussv-ont-pipeline
```

## Citation

If you use ConsensuSV-ONT, please cite:
ConsensuSV-ONT - a modern method for accurate structural variant calling
Antoni Pietryga, Mateusz Chiliński, Sachin Gadakh, Dariusz Plewczynski
bioRxiv 2024.07.26.605267; doi: https://doi.org/10.1101/2024.07.26.605267

## Preparation of your samples
Since the algorithm is easy-to-run, the preparation of the sample is minimal. The algorithm works with fastq.gz files, and a csv file containing all the sample information is needed. There is only one column, that is the path to the fastq.gz file.

Path |
-------------- |
/HG00733.fastq.gz |
/HG00514.fastq.gz |

Bear in mind that the column headers are provided only for the ease of the example, and should not be present in the csv file.

## Running the pipeline
To get into the docker image for working with your data, it's best to mount local directories to the container:

```bash
docker run -v /mnt:/mnt -p 8082:8082 -it antonipietryga/consensussv-ont-pipeline
```

Once in the container, we can run the pipeline using the following command:
```bash
nextflow pipeline.nf --input files.csv 
```

Remember to put your samples in the file samples.csv according to the [praparation of your samples](#preparation-of-your-samples).

All the parameters that can be used with the script are shown in the following table:

Parameter | Description
-------------- | ---------------
--input | File location of the csv file that described all the samples according to the (#preparation-of-your-samples).
--threads | Max number of threads per task.
--mem | Max memory per thread.
--outdir | Output dir of ConsensuSV-ONT.
--ref | Reference genome

## Output location

The location of the output depends on your working directory, provided as the parameter. In that directory, two folders will be created:
* alignments - a folder where alignment files for each sample are stored
* vcfs - where you will find folders for each of the sample, containing all VCF files from the individual SV calling tools and consensuSV-ONT

## CNN model architecture
```Model: "sequential"
_________________________________________________________________
Layer (type)                Output Shape              Param #   
=================================================================
 conv3d (Conv3D)             (None, 50, 50, 3, 4)      112       
                                                                 
 batch_normalization (BatchN  (None, 50, 50, 3, 4)     16        
 ormalization)                                                   
                                                                 
 max_pooling3d (MaxPooling3D  (None, 25, 25, 3, 4)     0         
 )                                                               
                                                                 
 conv3d_1 (Conv3D)           (None, 25, 25, 3, 8)      872       
                                                                 
 batch_normalization_1 (Batc  (None, 25, 25, 3, 8)     32        
 hNormalization)                                                 
                                                                 
 max_pooling3d_1 (MaxPooling  (None, 12, 12, 3, 8)     0         
 3D)                                                             
                                                                 
 conv3d_2 (Conv3D)           (None, 12, 12, 3, 16)     3472      
                                                                 
 batch_normalization_2 (Batc  (None, 12, 12, 3, 16)    64        
 hNormalization)                                                 
                                                                 
 max_pooling3d_2 (MaxPooling  (None, 6, 6, 3, 16)      0         
 3D)                                                             
                                                                 
 flatten (Flatten)           (None, 1728)              0         
                                                                 
 dense (Dense)               (None, 32)                55328     
                                                                 
 dense_1 (Dense)             (None, 1)                 33        
                                                                 
=================================================================
Total params: 59,929
Trainable params: 59,873
Non-trainable params: 56```
