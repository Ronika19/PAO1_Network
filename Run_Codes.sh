cd Acute_Burn_Wound
mkdir -p Output_Files
python3 DEG_NonDEG_Module_Pathway.py
cd ..
echo Acute_Burn_Wound Done

cd Chronic_Surgical_Wound
mkdir -p Output_Files
python3 DEG_NonDEG_Module_Pathway.py
cd ..
echo Chronic_Surgical_Wound Done

cd Antibiotic_Resistance
mkdir -p Output_Files
python3 DEG_NonDEG_Module_Pathway.py
cd ..
echo Antibiotic_Resistance Done

