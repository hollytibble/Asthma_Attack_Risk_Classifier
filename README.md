<!-- Grab your social icons from https://github.com/carlsednaoui/gitsocial -->

# Development and Validation of a Machine Learning Risk Prediction Model for Asthma Attacks in Adults in Primary Care


Welcome!  In this project,  we provide a robust generalizable methodology for mining large longitudinal electronic health care data towards predicting clinical outcomes, using asthma and asthma attacks as a testbed. We thoroughly investigate different statistical learning models with different configurations to optimize the developed models, and systematically validate out-of-sample performance to objectively report model generalizability. We envisage these results could be translated to clinical practice for asthma, and the reported methodology may be adapted to study other related chronic diseases on the basis of mining routinely collected longitudinal clinical data.

The protocol paper for this study is entitled ["Predicting asthma attacks in primary care: protocol for developing a machine learning-based prediction model"](https://bmjopen.bmj.com/content/9/7/e028375) (Tibble et al., 2019, 9:e028375 BMJ Open).

Thank you for visiting,

Holly Tibble [![alt text][1.2]][1]



## Paper Abstract
Introduction:
Primary care consultations provide an opportunity for patients and clinicians to assess asthma attack risk. Using a data-driven risk prediction tool with routinely collected health records may be an efficient way to aid promotion of effective self-management, and support clinical decision making. 

Methods:
Longitudinal Scottish primary care data for 21,250 asthma patients were used to predict the risk of asthma attacks in the following year.  A selection of machine learning algorithms (i.e., Naïve Bayes Classifier, Random Forests, and Extreme Gradient Boosting), hyperparameters, training data enrichment methods were explored, and validated in a random unseen data partition. 

Results:
Our final Extreme Gradient Boosting model achieved the best performance when no training data enrichment was applied. Around 1 in 3 (36.2%) predicted high-risk patients had an attack within one year of consultation, compared to approximately 1 in 16 in the predicted low-risk group (6.7%).  The model was well calibrated, with a calibration slope of 1.02 and an intercept of 0.004, and the Area under the Curve was 0.75.

Conclusion:
This model has the potential to increase the efficiency of routine asthma care by creating new personalized care pathways mapped to predicted risk of asthma attacks, such as priority ranking patients for scheduled consultations and interventions. Furthermore, it could be used to educate patients about their individual risk and risk factors, and promote healthier lifestyle changes, use of self-management plans, and early emergency care seeking following rapid symptom deterioration. 
port tool to enhance patient care. 


[1]: https://twitter.com/HollyTibble
[1.2]: http://i.imgur.com/wWzX9uB.png
