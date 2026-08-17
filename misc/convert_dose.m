dose = double(squeeze(dicomread("original_field_dose.dcm"))); 
info = dicominfo("original_field_dose.dcm"); 

dose = dose * info.DoseGridScaling; 
