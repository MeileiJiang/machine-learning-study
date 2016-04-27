######################################################################################
## NonlinearRegression.R
## Author: Meilei
## Date: April 2016
######################################################################################
library(ggplot2)
library(dplyr)
library(ISLR)
library(MASS)

Wage = Wage
glimpse(Wage)


# polynomial regression ---------------------------------------------------

ggplot(data = Wage, aes(x = age, y = wage)) +
  geom_point(col = "blue") +
  scale_y_continuous(limits = c(0, 350)) +
  geom_smooth(col = "red", formula = "y ~ poly(x, 4)") +
  labs(x = "Age", y = "Wage")


wage.poly <- lm(wage ~ poly(age, degree = 4), data = Wage)
summary(wage.poly)
wage.predict = predict(wage.poly)
plot(wage.predict~ Wage$age, type = "l")
