data <- readRDS("data2.rds")
genomemap <- readRDS("genomemap2.rds")
samplekey <- readRDS("samplekey2.rds")

library(randomForest)
library(datasets)
library(caret)
library(nnet)
library(pROC)
library(e1071)  # For SVM
library(class)  # For KNN
library(gbm)    # For GBM
library(MASS)   # For LDA
### Random Forest Classification ##############################################
set.seed(123)
data_variable <- data[order(apply(data, 1, sd),
                      decreasing = TRUE),][1:1000, ]

ind <- sample(2, ncol(data_variable), replace = TRUE, prob = c(0.8, 0.2))
train <- data_variable[,ind==1]
test <- data_variable[,ind==2]

train_transposed <- as.data.frame(t(train))
train_donors <- row.names(train_transposed)
train_transposed$smoking <- samplekey[train_donors,]$smoking
train_transposed$smoking <- as.factor(train_transposed$smoking)
model_rf <- randomForest(smoking ~ ., data = train_transposed,ntree = 1000, proximity=TRUE)
print(model_rf)


test_transposed <- as.data.frame(t(test))
test_donors <- row.names(test_transposed)
test_transposed$smoking <- samplekey[test_donors,]$smoking
test_transposed$smoking <- as.factor(test_transposed$smoking)
pr <- predict(model_rf, test_transposed)
confusionMatrix(pr, test_transposed$ smoking)

plot(model_rf$err.rate[,1], type='l', col='blue', ylim=c(0,0.5), 
     ylab='Error Rate', main='Error Rate by Number of Trees', xlab='Number of Trees')
lines(model_rf$err.rate[,2], type='l', col='red')
legend('topright', legend=c('Training', 'Out-of-Bag'), col=c('blue', 'red'), lty=1, y = 0)




### Validation #################################################################
train_control <- trainControl(method = "cv", number = 10)
classification_model <- train(smoking ~ ., data = train_transposed,
                              method = "rf", trControl = train_control)
print(classification_model)
pr_class <- predict(classification_model, newdata = test_transposed)
confusionMatrix(pr_class, test_transposed$smoking)

####################
# Function to evaluate classifiers
evaluate_classifier <- function(model, test_data, test_labels) {
  predictions <- predict(model, newdata = test_data)
  conf_matrix <- confusionMatrix(predictions, test_labels)
  return(list(predictions = predictions, confusion_matrix = conf_matrix))
}

# Train and evaluate Random Forest
model_rf <- train(smoking ~ ., data = train_transposed,
                method = "rf", trControl = train_control)
rf_results <- evaluate_classifier(model_rf,
                  test_transposed, test_transposed$smoking)

# Train and evaluate SVM
model_svm <- train(smoking ~ ., data = train_transposed,
                  method = "svmRadial", trControl = train_control)
svm_results <- evaluate_classifier(model_svm,
                  test_transposed, test_transposed$smoking)

# Train and evaluate KNN
model_knn <- train(smoking ~ ., data = train_transposed, method = "knn",
                   trControl = train_control, tuneLength = 10)
knn_results <- evaluate_classifier(model_knn,
                  test_transposed, test_transposed$smoking)

# Train and evaluate GBM
model_gbm <- train(smoking ~ ., data = train_transposed, method = "gbm",
                  trControl = train_control, verbose = FALSE)

model_gbm_predict <- predict(model_gbm, newdata = test_transposed)
confusionMatrix(model_gbm_predict, test_transposed$smoking)
gbm_results <- evaluate_classifier(model_gbm,
                  test_transposed, test_transposed$smoking)

# Train and evaluate Logistic Regression
model_glm <- train(smoking ~ ., data = train_transposed, method = "glm",
                   family = binomial, trControl = train_control)
glm_results <- evaluate_classifier(model_glm, test_transposed,
                    test_transposed$smoking)

# Train and evaluate Linear Discriminant Analysis (LDA)
model_lda <- train(smoking ~ ., data = train_transposed,
                   method = "lda", trControl = train_control)
lda_results <- evaluate_classifier(model_lda,
                    test_transposed, test_transposed$smoking)

results <- data.frame(
  Model = c("Random Forest", "SVM", "KNN", "GBM", "Logistic Regression", "LDA"),
  Accuracy = c(rf_results$confusion_matrix$overall['Accuracy'], 
               svm_results$confusion_matrix$overall['Accuracy'], 
               knn_results$confusion_matrix$overall['Accuracy'], 
               gbm_results$confusion_matrix$overall['Accuracy'], 
               glm_results$confusion_matrix$overall['Accuracy'], 
               lda_results$confusion_matrix$overall['Accuracy']),
  Kappa = c(rf_results$confusion_matrix$overall['Kappa'], 
            svm_results$confusion_matrix$overall['Kappa'], 
            knn_results$confusion_matrix$overall['Kappa'], 
            gbm_results$confusion_matrix$overall['Kappa'], 
            glm_results$confusion_matrix$overall['Kappa'], 
            lda_results$confusion_matrix$overall['Kappa'])
)

print(results)


### Principine komponente ######################################################
pc <- prcomp(train_transposed[,-1001],
             center = TRUE,
             scale. = TRUE)
trg <- predict(pc, train_transposed)
trg <- data.frame(trg, train_transposed[1001])
tst <- predict(pc, test_transposed)
tst <- data.frame(tst, test_transposed[1001])

pca_model <- multinom(smoking~., data = trg)

print(pca_model)
p <- predict(pca_model, trg)
confusionMatrix(p, trg$smoking)
p1 <- predict(pca_model, tst)
confusionMatrix(p1, tst$smoking)
