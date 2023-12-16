


#Rscript ./process_beam_bcr.R --rawdir "./23_303_3/" --save_name "test" --lanes 1,2
#Rscript ./process_beam_bcr4rawdata.R --rawdir "./RAW/" --save_name "raw_test" --lanes 1,2


seur <- dsbNormalizeSeurat("./raw_test_20231216012515/raw_test_20231216012515.RDS", "./test_20231216012454/test_20231216012454.RDS")
seurat_obj <-readRDS("./test_20231216012454/test_20231216012454.RDS")


