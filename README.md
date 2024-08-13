# Advanced-topics-in-physiological-signal-processing
This project aims to replicate the results of the study titled "Paroxysmal Slow Wave Events Predict Epilepsy Following a First Seizure," published in 'Epilepsia'. The original study developed a logistic regression model based on four key features: occurrence per minute, mean madian power frequency (MPF), mean duration, and the average number of channels detecting PSWEs. For replication, I conducted a logistic regression analysis with a cohort of 338 patients, comprising 184 with epilepsy and 154 without. Expanding upon the original study, I included six additional features and used a forward feature selection (FFS) algorithm to identify the most informative features, aiming to reduce overfitting and improve model generalization.

The FFS algorithm identified four of the ten features as optimal: mean MPF, mean relative power in delta, mean relative power in beta, and mean event duration. This model demonstrated strong performance in distinguishing between epilepsy and non-epilepsy cases when evaluated with validation data. In contrast, a logistic regression model trained with only the original four features showed inferior area under the curve (AUC), sensitivity, and specificity. The model incorporating FFS also had a higher maximum odds ratio (30.9934) compared to the odds ratio of 8.79 in the replication of the article, indicating a stronger predictive power for epilepsy-related features.

This higher odds ratio suggests that the selected features or feature selection methods in my model are more effective at predicting epilepsy. The significant differences observed highlight the importance of feature selection in model performance. However, the variability in results across different runs suggests that model performance may fluctuate due to factors such as data variability, outliers, or noise. To further understand and improve the model's reliability, additional analysis, including increasing the dataset size and exploring alternative modeling approaches, may be necessary. Additional results are provided in the appendices.

![image](https://github.com/user-attachments/assets/5d1c4236-af84-4e26-83d0-6c89ebadc9b5)
Figure 1:  FFS results - The AUC as function of the number of chosen features

![image](https://github.com/user-attachments/assets/b23ca46f-5814-48f4-a53f-98f9752c2032)
Figure 2:  The ROC of the chosen model on the test data. 
In red- the threshold that chosen using Youden J statistic.


