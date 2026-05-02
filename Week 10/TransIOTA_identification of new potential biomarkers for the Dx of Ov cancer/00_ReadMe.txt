-------------------------------------------
GENERAL INFORMATION
-------------------------------------------

G01. Names of file(s) or dataset(s) that this README file describes

This file describes the datasets "TransIOTA - raw data.csv" and "TransIOTA - Imputations.xlsx".
"TransIOTA - raw data.csv" contains all the data collected during the TransIOTA-study.
"TransIOTA - Imputations.xlsx" contains the imputed values for variables with missing observations (protein measurements, immune cells measurements and cell-free DNA scores). Imputations for different analyses are shown in seperate sheets.
The datasets are accompanied with a data dictionary 'TransIOTA - data dictionary.xlsx' and a masterdocument 'Masterdocument TransIOTA - RDR.xlsx'.
The data dictionary gives additional information about the variables.
The masterdocument gives an overview of the patients used in each analyses.


G02. Date of creation/last update of the README file

Created on 2023-August-28
Last updated on 2023-October-9


G03. Name and contact information of Principal Investigator

Prof. An Coosemans (an.coosemans@kuleuven.be)


G04. ORCID of Principal Investigator

0000-0002-7321-4339


G05. Institution of Principal Investigator

KU Leuven, BELGIUM


G06. Contact of other person at KU Leuven that has access to the dataset

Prof. Ben Van Calster (ben.vancalster@kuleuven.be)

Prof. Dirk Timmerman (dirk.timmerman@uzleuven.be)


G07. Description of the dataset

The dataset contains clinical information, ultrasound findings and blood measurements (immune cells, proteins and cell-free DNA) of patients with ovarian masses that underwent surgery. 


G08. Keywords (author defined)

TransIOTA, biomarkers, cfDNA, genome-wide z-score, nucleosome score, prediction models, ROMA, ADNEX, validation


G09. Thesaurus or controlled vocabulary keywords

Ovarian neoplasms (Mesh)


G10. Thesaurus or controlled vocabulary used in this README

/


G11. Language(s) used in the dataset

English


G12. Other involved researchers

Jolien Ceusters (ORCID: 0000-0001-5422-5668) 


------------------------------------------
PROJECT INFORMATION
------------------------------------------

P01. Project information

This data is collected for the TRANS-IOTA (translational IOTA) research project. The aim of TRANS-IOTA is to search for new biomarkers in ovarian cancer patients that can serve as new diagnostic tool on itself or can be integrated in the IOTA models. 


P02. Project abstract

The TRANS-IOTA study looked at different biomarkers: protein measurements, immune cells and cell-free DNA scores. We investigated the potential of these biomarkers, on itself or in combination with ultrasound features, to improve the preoperative differential diagnosis between benign and malignant ovarian masses.
This data is also used to validate the performance of two prediction models: the Assessment of Different NEoplasias in the adneXa (ADNEX) model and the Risk of Ovarian Malignancy Algorithm (ROMA).


P03. Project funder: Name of funder, type of grant, grant number

Kom op Tegen Kanker (Stand up to Cancer): 2016/10728/2603
 

------------------------------------------
FILE OVERVIEW
------------------------------------------

F01. Number of files described by the README-file

4


F02. List with names of files, description, date of creation of file

TransIOTA - raw data.csv   -   Clinical information, ultrasound findings and blood measurements   -   9/10/2023
TransIOTA - Imputations.xlsx   -   Imputed values for missing measurements of the biomarkers   -   9/10/2023
TransIOTA - data dictionary.xlsx   -   Description of the variables    -   9/10/2023
Masterdocument TransIOTA - RDR.xlsx   -   Overview of patients used in the different papers   -   9/10/2023


F03. File formats

.csv, .xlsx


F04. Software used to generate the data

Clinical and ultrasound data is collected with Clinical Data Miner, a web-based software with an electronic Case Report Form (eCRF).
Blood measurements (cfDNA scores, proteins and immune cells) were collected in Excel-files.
Imputations were created with R.
The data dictionary and masterdocument are created with Excel.


F05. Software necessary to open the file

Excel
 

F06. Relationship between the files

'TransIOTA - raw data.csv' and 'TransIOTA - Imputations.xlsx' can be linked with each other based on the variable 'Patient.ID'.


F07. Which version of the dataset is this? Date of this version?

First version


F08. Information about the dataset versions and reason for updates

/


F09. Naming conventions for file names

/


--------------------------------------------
STORAGE INFORMATION
--------------------------------------------

S01. Where are the data stored?

- KU Leuven RDR (this repository)


S02. Links to other available locations of the dataset (e.g. repository)

/


--------------------------------------------
METHODOLOGICAL INFORMATION
--------------------------------------------
M01. Date (beginning-end) and place of data collection

June 2015 - December 2019: University Hospitals Leuven (Leuven, Belgium), Università Cattolica del Sacro Cuore (Rome, Italy), General Faculty Hospital of the Charles University (Prague, Czech Republic), Queen Charlotte’s and Chelsea Hospital, Imperial College (London, UK), Ziekenhuis Oost-Limburg (Genk, Belgium), and Istituto Nazionale dei Tumore (Milan, Italy).


M02. Aim for which the data were collected

Research


M03. Data collecting method

Clinical Data Miner (CDM)


M04. Information about data processing methods

The statistical analyses were done in R. Details about the performed analyses can be found in the manuscripts.


M05. Information about the instrument, calibration

/


M06. Quality assurance procedures

/


M07. Information about limitations of the dataset, information that ensures correct interpretation of the dataset

The dataset contains only patients with an ovarian mass that underwent surgery within 120 days after the recruitment ultrasound scan.


M08. People involved in the creation or processing of the dataset



--------------------------------------------
DATA ACCESS AND SHARING
--------------------------------------------

A01. Recommended citation for the dataset

/


A02. License information, restrictions on use

The dataset is not publicly available. However, the dataset may be obtained following permission of the contact persons (Prof. An Coosemans and Prof. Dirk Timmerman)


A03. Confidentiality information

/


----------------------------------------------------------
DATA SPECIFIC INFORMATION (ABOUT THE DATA THEMSELVES)
----------------------------------------------------------

D01. Full names and definitions for columns and rows

See data dictionary


D02. Explanation of abbreviations

See data dictionary


D03. Units of measurement

Maximum diameter of the lesion - mm
Maximum diameter of the solid component - mm
Serum CA125 - kU/L
Serum HE4 - pmol/L
More details in the data dictionary.


D04. Symbols for missing data

OOR = out of range: This are protein measurements that were outside of the detection limit.


-----------------------------------------------------
RELATIONSHIPS
-----------------------------------------------------

R01. Publications based on this dataset

- Coosemans, A., Baert, T., Ceusters, J., Busschaert, P., Landolfo, C., Verschuere, T., ... & Timmerman, D. (2019). Myeloid-derived suppressor cells at diagnosis may discriminate between benign and malignant ovarian tumors. International Journal of Gynecologic Cancer, 29(9).
  Included cases are shown in column I (Immune cells) of the masterdocument
- Landolfo, C., Achten, E. T. L., Ceusters, J., Baert, T., Froyman, W., Heremans, R., ... & Coosemans, A. (2020). Assessment of protein biomarkers for preoperative differential diagnosis between benign and malignant ovarian tumors. Gynecologic Oncology, 159(3), 811-819.
  Included cases are shown in column J (Proteins (interim)) of the masterdocument
- Vanderstichele, A., Busschaert, P., Landolfo, C., Olbrecht, S., Coosemans, A., Froyman, W., ... & Vergote, I. (2022). Nucleosome footprinting in plasma cell-free DNA for the pre-surgical diagnosis of ovarian cancer. NPJ Genomic Medicine, 7(1), 30.
  Included cases are shown in column K (Nucleosome) of the masterdocument
- "Comparison of the ADNEX and ROMA risk prediction models for the diagnosis of ovarian cancer: a multicentre external validation in patients who underwent surgery" - Submitted for review
  Included cases are shown in column H (ROMA vs ADNEX) of the masterdocument
- "Added value of proteins to clinical and ultrasound information in predicting the risk of malignancy in ovarian tumors" - In preparation
  Included cases are shown in column G (Proteins) of the masterdocument
- "Added value of plasma cell-free DNA to the ADNEX model for predicting the risk of ovarian cancer" - In preparation
  Included cases are shown in column F (cfDNA) of the masterdocument


R02. This dataset derives from (other dataset)

/


R03. This dataset is related to (documents, dataset)

/


R04. References of publications used to create the datasets

/



