# ConsensusSV-ONT-pipeline

# CNN model architecture
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
