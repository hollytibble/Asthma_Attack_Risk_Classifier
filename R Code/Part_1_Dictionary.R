#####################################################################
###  Read Code Dictionary
#####################################################################

### Asthma Diagnosis
codes_asthma_diagnosis <- c("173A.", "H3120", "H33..", "H330.", "H3300", "H3301", "H330z", "H331.",
                            "H3310", "H3311", "H331z", "H332.", "H334.", "H335.", "H33z.", "H33z0", 
                            "H33z1", "H33z2", "H33zz", "H3B..", "663..", "6632.", "6636.", "6637.", 
                            "663a.", "663B.", "663d.", "663e.", "663e0", "663e1", "663F.", "663f.", 
                            "663G.", "663g.", "663g0", "663g1", "663g2", "663g3", "663H.", "663h.", 
                            "663I.", "663J.", "663j.", "663L.", "663M.", "663m.", "663N.", "663n.",
                            "663N0", "663N1", "663N2", "663O.", "663O0", "663P.", "663p.", "663Q.", 
                            "663q.", "663R.", "663r.", "663S.", "663s.", "663T.", "663t.", "663U.", 
                            "663u.", "663V.", "663v.", "663V0", "663V1", "663V2", "663V3", "663W.",
                            "663w.", "663X.", "663x.", "663Y.", "663y.", "663Z.", "663z.")


### Asthma Related Codes
codes_asthma_related<-c("173c.", "173d.", "178..", "1O2..", "388t.", "66Y0.", "66Y1.",
                        "66Y2.", "66Y3.", "66Y4.", "66Y5.", "66Y6.", "66Y7.", "66Y8.", 
                        "66Y9.", "66YA.", "66Ya.", "66Yb.", "66YC.", "66Yc.", "66YE.", 
                        "66YF.", "66YG.", "66YJ.", "66YK.", "66Ym.", "66YN.", "66YO.",
                        "66YP.", "66YQ.", "66YR.", "66YV.", "66YW.", "66YX.", "66YY.",
                        "66YZ.", "679J.", "8B3j.", "8CE2.", "8CR0.", "8HTT.", "9N1d.", 
                        "9NI8.", "9OJ..", "9OJ1.", "9OJA.", "H334.", "H33z.", "H33z2",
                        "H33zz")


### Asthma Attack
codes_asthma_attack<-c("H3301", "H3311", "H33z0", "H33z1", "H333.",
                       "663d.", "663m.", "8H2P.", "663y.")

### COPD
codes_COPD<-c("H3...","H36..","H37..","H38..","H39..",
              "H3A..", "H3y1.", "H3y..", "H3z..", 
              "H31..", "H310.", "H3100", "H310z",
              "H311.", "H3110", "H3111", "H311z", "H312.",
              "H3121", "H3123","H312z", "H313.",
              "H31y.", "H31y1","H31yz", "H31z.", "H32..",
              "H320.", "H3200", "H3201","H3202", "H3203", 
              "H320z", "H321.", "H322.", "H32y.", "H32y0", 
              "H32y1", "H32y2", "H32yz", "H32z.")

### Nasal Polpys
codes_nasal_polyps<-c("H11..")


### Smoking
codes_smoking_never <- c("1371.")
codes_smoking_former<-c("1377.", "1378.", "1379.", "137A.", "137B.",
                        "137F.", "137i.", "137j.", "137K.", "137K0",
                        "137L.", "137l.", "137N.", "137O.", "137S.",
                        "137T.")
codes_smoking_current<-c("1372.", "1373.", "1374.", "1375.", "1376.",
                         "137a.", "137b.", "137c.", "137C.", "137d.",
                         "137D.", "137e.", "137f.", "137G.", "137h.",
                         "137H.", "137J.", "137M.", "137m.", "137P.",
                         "137Q.", "137R.", "137V.")
codes_smoking_current_conditional<-c("137..", "137E.", "137g.", "137X.", "137Y.",
                                     "137Z.")
codes_smoking_all<-c(codes_smoking_never,
                     codes_smoking_former,
                     codes_smoking_current, 
                     codes_smoking_current_conditional)


### BMI
codes_BMI_val<-c("22K..") 
codes_BMI_under<-c("22K3.","22K6.")
codes_BMI_normal<-c("22K1.","22K8.")
codes_BMI_overweight<-c("22K2.","22K4.") # BMI high coded here as overweight rather than obese
codes_BMI_obese<-c("22K5.","22K7.", "22KC.","22KD.", "22KE.")
codes_BMI_cat<-c(codes_BMI_under,codes_BMI_normal,codes_BMI_overweight,codes_BMI_obese)
codes_height<-"229.."
codes_weight<-"22A.."

### Rhinitis
codes_rhinitis<-c("H17..", "H170.", "H171.", "H1710", "H172.",
                  "H17z.", "H18..", "Hyu21")

### Eczema / Dermitits
codes_eczema<-c("M11..", "M111.", "M112.","M113.", "M114.", "M11z.",
                "M12z0", "M12z1")

### GERD
codes_GERD<-c("J101.","J10y4","J10y6","J1011","J1016",
              "J101z", "J1025","J1020","1957.")

### Blood Eosinophil Counts
codes_eosinophils<-"42K.."

### Anxiety and Depression
codes_anx_dep<-c("8G94.", "E2...", "E20..", "E200.", "E2000", "E2001", "E2002", "E2003", 
                 "E2004", "E2005", "E200z", "E201.", "E2010", "E2011", "E2012", "E2013", 
                 "E2014", "E2015", "E2016", "E2017", "E2018", "E2019", "E201A", "E201B", 
                 "E201C", "E201z", "E202.", "E2020", "E2021", "E2022", "E2023", "E2024", 
                 "E2025", "E2026", "E2027", "E2028", "E2029", "E202A", "E202B", "E202C", 
                 "E202D", "E202E", "E202z", "E203.", "E2030", "E2031", "E203z", "E205.", 
                 "E206.", "E207.", "E20y.", "E20y0", "E20y1", "E20y2", "E20y3", "E20yz", 
                 "E20z.", "E21..", "E210.", "E211.", "E2110", "E2111", "E2112", "E2113", 
                 "E211z", "E26..", "E260.", "E2600", "E2601", "E260z", "E261.", "E2610", 
                 "E2611", "E2612", "E2613", "E2614", "E2615", "E261z", "E262.", "E2620", 
                 "E2621", "E2622", "E2623", "E262z", "E263.", "E2630", "E263z", "E264.", 
                 "E2640", "E2642", "E2643", "E2644", "E2645", "E264z", "E265.", "E2650", 
                 "E2651", "E2652", "E2653", "E265z", "E266.", "E267.", "E26y.", "E26y0", 
                 "E26yz", "E26z.", "E278.", "E2780", "E2781", "E2782", "E278z", "E28..", 
                 "E280.", "E281.", "E282.", "E283.", "E2830", "E2831", "E283z", "E284.", 
                 "E28z.", "E29..", "E2900", "E292.", "E2920", "E2921", "E2922", "E2923", 
                 "E2924", "E2925", "E292y", "E292z", "E293.", "E2930", "E2931", "E2932", 
                 "E293z", "E294.", "E29y.", "E29y0", "E29y1", "E29y2", "E29y3", "E29y4", 
                 "E29y5", "E29yz", "E29z.", "Eu4..", "Eu40.", "Eu400", "Eu401", "Eu402", 
                 "Eu403", "Eu40y", "Eu40z", "Eu41.", "Eu410", "Eu411", "Eu412", "Eu413", 
                 "Eu41y", "Eu41z", "Eu42.", "Eu420", "Eu421", "Eu422", "Eu42y", "Eu42z", 
                 "Eu43.", "Eu430", "Eu431", "Eu432", "Eu43y", "Eu43z", "Eu44.", "Eu440", 
                 "Eu441", "Eu442", "Eu443", "Eu444", "Eu445", "Eu446", "Eu447", "Eu44y", 
                 "Eu44z", "Eu45.", "Eu450", "Eu451", "Eu452", "Eu453", "Eu454", "Eu455", 
                 "Eu45y", "Eu45z", "Eu46.", "Eu460", "Eu461", "Eu46y", "Eu46z", "ZN114", 
                 "ZS7C7", "1B17.", "62T1.", "6G00.", "8CAa.", "9H90.", "9H91.", "9H92.", 
                 "E03y2", "E03y3", "E11..", "E112.", "E1120", "E1121", "E1122", "E1123", 
                 "E1124", "E1125", "E1126", "E112z", "E113.", "E1130", "E1131", "E1132", 
                 "E1133", "E1134", "E1135", "E1136", "E1137", "E113z", "E118.", "E11y2", 
                 "E11y3", "E11yz", "E11z.", "E11z0", "E11z1", "E11z2", "E11zz", "E135.", 
                 "E204.", "E290.", "E290z", "E291.", "E2B..", "E2B0.", "E2B1.", "Eu3..", 
                 "Eu32.", "Eu320", "Eu321", "Eu322", "Eu324", "Eu32y", "Eu32z", "Eu33.", 
                 "Eu330", "Eu331", "Eu332", "Eu334", "Eu33y", "Eu33z", "Eu34.", "Eu340", 
                 "Eu341", "Eu34y", "Eu34z", "Eu3y.", "Eu3y0", "Eu3y1", "Eu3yy", "Eu3z.")

### Charlson Comorbidity Index
codes_charlson_pulmonary<- c("1761.", "1780.", "14B4.", "H30..", "H300.", "H30z.", 
                             "H3101", "H34..", "H340.", "H341.", "H34z.", "H35..", "H350.", 
                             "H351.", "H352.", "H3520", "H3521", "H352z", "H353.", "H354.", 
                             "H355.", "H356.", "H35y.", "H35y3", "H35y5", "H35y6", "H35y7", 
                             "H35yz", "H35z.", "H35z1", "H35zz", "H40..", "H41..", 
                             "H410.", "H41z.", "H42..", "H420.", "H423.", "H42z.", "H43..", 
                             "H430.", "H431.", "H432.", "H434.", "H435.", "H43z.", "H440.", 
                             "H441.", "H442.", "H45..", "H460.", "H460z", "H4640", "H4641", 
                             "H4642", "H47y0", "H4y10", "H4z..", "H57y.", "H57yz", "H581.", 
                             "H582.", "Hyu30", "Hyu40", "Hyu41", "Hyu43", "SK07.")
  
### Peak Flow
codes_peak_flow<-c("339A.", "339c.")

### Respiratory Infection
codes_RI<-c("H0...","H05..","H05z.","H06..","H06z.","H07..","H0y..","H04..","H05y.","H0z..",
            "H060.","H0603","H0604","H0605","H0609","H060D","H060X","H060w","H060z","H061.",
            "H0612","H0615","H0616","H061z","H062.","H06z0","H06z1","H06z2", "H2...","H20..",
            "H201.","H20y.","H20y0","H20z.","H21..","H22..","H220.","H222.","H223.","H224.",
            "H22y.","H22y2","H22yz","H22z.","H23..","H231.","H23z.","H24.","H24y.","H24y2",
            "H24yz","H24z.","H25..","H26..","H260.","H2600","H261.","H262.","H263.","H27..",
            "H270.","H2701","H271.","H2710","H27z.","H28..","H2A..","H2B..","H2C..","H2y..",
            "H2z..","H300.", "H301.", "H30z.", "H5400", "H5401","Hyu0.","Hyu08","Hyu0A",
            "Hyu0b","Hyu0H","Hyu1.", "Hyu10", "Hyu11","G5203")
            

#####################################################################
###  Prescribing Keywords Dictionary
#####################################################################

### Generic Drugs
med_SABA<-c("SALBUTAMOL")
med_LABA<-c("BAMBUTEROL", "FORMOTEROL", "SALMETEROL", "TERBUTALINE", "TIOTROPIUM")
med_LAMA<-c("IPRATROPIUM")
med_Theophylline<-c("THEOPHYLLINE", "AMINOPHYLLINE")
med_ICS<-c("BECLOMETASONE", "BUDESONIDE", "FLUTICASONE", 
           "MOMETASONE", "CICLESONIDE")
med_LTRA<-c("MONTELUKAST", "ZAFIRLUKAST", "CROMOLYN", "NEDOCROMIL")
med_Steroid<-c("PREDNISOLONE","PREDNISOLONE","DEXAMETHASONE")
med_MAb<-c("OMALIZUMAB", "MEPOLIZUMAB", "BENRALIZUMAB", "RESLIZUMAB")
### List of all medicine keywords
med_list<-substr(ls()[which(substr(ls(),1,3)=="med")],5,30)

### Brand named drugs
brand_SALBUTAMOL<-c("SALAMOL", "VENTOLIN", "AIROMIR", "SALBULIN",
                    "AIRSALB", "VENTMAX", "ASMASAL", 
                    "PULVINAL SALBUTAMOL","ALBUTEROL")
brand_BAMBUTEROL<-c("BAMBEC")
brand_FORMOTEROL<-c("ATIMOS", "FORADIL", "OXIS", "FOSTAIR",
                    "SYMBICORT", "DUORESP SPIROMAX", "FOBUMIX", "FLUTIFORM")
brand_SALMETEROL<-c("NEOVENT", "SEREVENT", "SOLTEL","SERETIDE", "AIRFLUSAL", 
                    "SIRDUPLA", "SEREFLO", "ALOFLUTE", "COMBISAL", 
                    "FUSACOMB", "STALPEX")
brand_TERBUTALINE<-c("BRICANYL")
brand_TIOTROPIUM<-c("SPIRIVA RESPIMAT")
brand_IPRATROPIUM<-c("ATROVENT", "INHALVENT", "IPRAVENT", "RESPONTIN", "IPRAMOL", 
                     "COMBIVENT")
brand_THEOPHYLLINE<-c("NUELIN", "SLOPHYLLIN", "SLO-PHYLLIN","UNIPHYLLIN")
brand_AMINOPHYLLINE<-c("PHYLLOCONTIN")
brand_BECLOMETASONE<-c("BECLOMETHASONE","CLENIL", "QVAR", "KELHALE", "SOPROBEC", 
                       "BECODISKS", "PULVINAL BECLOMETASONE", "ASMABEC","FOSTAIR")
brand_BUDESONIDE<-c("BUDELIN", "PULMICORT","SYMBICORT", "DUORESP SPIROMAX", "FOBUMIX")
brand_MOMETASONE<-c("ASMANEX", "TWISTHALER")
brand_FLUTICASONE<-c("FLIXOTIDE","FLUTIFORM", "SERETIDE", "AIRFLUSAL", 
                     "SIRDUPLA", "SEREFLO", "ALOFLUTE", "COMBISAL", 
                     "FUSACOMB", "STALPEX","RELVAR ELLIPTA")
brand_CICLESONIDE<-c("ALVESCO")
brand_MONTELUKAST<-c("SINGULAIR")
brand_ZAFIRLUKAST<-c("ACCOLATE")
brand_CROMOLYN<-c("CROMOGLICATE", "INTAL")
brand_NEDOCROMIL<-c("TILADE")
brand_PREDNISOLONE<-c("DELATCORTRIL","PEVANTI","DILACORT","DELTASTAB")
brand_DEXAMETHASONE<-c()
brand_PREDNISONE<-c()
brand_OMALIZUMAB<-c("XOLAIR")
brand_MEPOLIZUMAB<-c("NUCALA")
brand_BENRALIZUMAB<-c("FASENRA")
brand_RESLIZUMAB<-c("CINQAERO")

### List of all medicine keywords
brand_list<-substr(ls()[which(substr(ls(),1,5)=="brand")],7,30)

exclusion_brands<-c("NASONEX","FLIXONASE","ANORO ELLIPTA","SUMATRIPTAN","AVAMYS",
                    "RHINOCORT","NASOBEC","NASOFAN","RYNACROM","PIRINASE","SPIOLTO",
                    "DYMISTA","POLLENASE","VIVIDRIN","DUAKLIR","SEEBRI", "ULTIBRO", 
                    "PRED FORTE", "TRELEGY", "TRIMBOW", "BRALTUS", "RINATEC", 
                    "ENTOCORT", "BENACORT", "AIRCORT", "BUDEFLAM", "BUDENOFALK",
                    "CORTIMENT", "JORVEZA", "AZELASTINE", "CUTIVATE", "ELOCON", 
                    "NALCROM", "CATACROM", "ASPIRE", "OPTICROM", "OPTREX", "BECONASE",
                    "MURINE", "ACLIDINIUM", "GENUAIR", "OLADATEROL", "YANIMO")

exclusion_keywords<-c("NASAL","NOSE","NOSTRIL","NASULE", "HAYFEVER",
                          "EYE","EAR","DROP","TONGUE", "FOAM","ENEMA", 
                          "RECTAL","SUPPOSITOR", "CREAM",
                          "OINTMENT", "ULCER","SKIN","PATCH","APPLY")

steroid_ind<-c("FOR MS","VASCULITIS","MYSTHAENIA","TRANSPLANT","COLITIS",
               "TISSUE DISEASE","CROHNS","PMR","RENAL","MYOSITIS",
               "MUSCLE","GOUT","ARTERITIS","MOBILITY","URTICARIA","COPD",
               "CROUP","IMMUN","HEADACHE","HEPATITIS","SKIN","RASH",
               "MYALGIA","RHEUM","ARTHRITIS","PAIN","DERMATOLOG",
               "SWELLING","EBLOW","SHOULDER","LUPUS","SCLEROSIS","ECZEMA",
               "HIVES","BOWEL","SCIATICA","HAYFEVER","RHINITIS","JOINT","NERVE")

### formulation keywords
formulation_keywords_sol<-c("SACHET", "RESPULE", "NEB", "VIAL", "AMPOULE")

### Doses
doses_mcg<-c("10000","5000","4000","2000","1000","500","400",
             "320","250","200","184","160","125","100","92","80","65","50")
doses_mg<-c("0.5","20","10","5","4","2","1")

### Pack Size Doses
pack_size_doses<-c(200,120,112,100,60,56,50,40,30,28,24,20,14,5)

### Dose Frequency Keywords
dose_freq_one<-c("ONCE","O-D","O\\.D")
dose_freq_two<-c("TWICE","TWO TIMES","2 TIMES","TD","TID","BID","BD",
                   "B-D","B\\.D")
dose_freq_four<-c("QID","FOUR TIMES","4 TIMES")
dose_freq_daily<-c("DAILY","EVERY DAY","EACH DAY",
                   "MANE","NOCTE","MORN","NIGHT","EVE","BEDTIME",
                   "A\\.M","P\\.M","AM","PM")

### Dose Quantity Keywords
dose_quant<-c("PUF","DOSE","CLICK","BLISTER", "TAB", "SACHET", "NEB", "RESP",
              "VIAL","CAP","INHALATION","AMPOULE","DOSE","ACTUATION","TWIST")


#####################################################################
###  Save Workspace
#####################################################################

save.image("../Data/Dictionary.RData")
