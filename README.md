<!-- Grab your social icons from https://github.com/carlsednaoui/gitsocial -->

# Short-Term Adult Asthma Attack Prediction using Electronic Health Record Data in the Primary Care Setting

Welcome!  In this project,  we provide a robust generalizable methodology for mining large longitudinal electronic health care data towards predicting clinical outcomes, using asthma and asthma attacks as a testbed. We thoroughly investigate different statistical learning models with different configurations to optimize the developed models, and systematically validate out-of-sample performance to objectively report model generalizability. We envisage these results could be translated to clinical practice for asthma, and the reported methodology may be adapted to study other related chronic diseases on the basis of mining routinely collected longitudinal clinical data.

The protocol paper for this study is entitled ["Predicting asthma attacks in primary care: protocol for developing a machine learning-based prediction model"](https://bmjopen.bmj.com/content/9/7/e028375) (Tibble et al., 2019, 9:e028375 BMJ Open).

The code, written in the R programming language, and clinical code lists for this project will be made availble in this repository once the paper has been published. 

Thank you for visiting,

Holly Tibble [![alt text][1.2]][1]



## Paper Abstract
Primary care consultations provide an opportunity for patients and clinicians to assess asthma attack risk.  Accurate prediction of that risk can prompt timely primary care intervention, increase regularity of primary care visits, promote risk-reducing lifestyle choices, and encourage patients to seek emergency care following symptom deterioration.  In this paper, we utilise the wealth of data recorded in electronic health records to develop and test a model for predicting asthma attacks within twelve weeks. 

The model features were patient demographics, smoking status, obesity, asthma treatment regimen and adherence, peak expiratory flow, history of asthma attacks, respiratory infections, blood eosinophil counts, and comorbidities.  The final analysis dataset comprised 40 features and 746,557 samples (asthma or respiratory infections primary care encounters) for 20,846 individuals, over eight years.

A repeated, random, split-sample approach was used to train and test multiple binary classification models, comparing imbalanced data handling approaches and statistical learning algorithms.  A 10% partition of the data samples was held-out for final model testing, with the remaining 90% used for model training, internal selection, and validation.  

In the unseen held-out data partition, our optimised random forest model had a balanced accuracy of 75.9% (sensitivity 52.3%, and specificity 99.6%).  Those with previous asthma attacks have higher sensitivity (66.9 versus 40.8%), with modest differences in positive predictive value (76.7% versus 74.5%).    

Our model can predict asthma attack incidence following an asthma or respiratory infection primary care consultation with high performance.  This can be developed into a clinical decision support tool to enhance patient care. 


[1]: https://twitter.com/HollyTibble
[1.2]: http://i.imgur.com/wWzX9uB.png
