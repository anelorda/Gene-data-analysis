---
title: "Data Analysis Project"
author: "Anel Ordabayeva"
date: "`r format(Sys.time(), '%Y-%m-%d')`"
output:
  html_document:
    toc: true
    toc_depth: 2
---

# Data Analysis test

The data frame “gene_tr” contains information about \~16,000 human genes:

• gene_id (character): Entrez gene id

• symbol (character): gene symbol

• cds_length (numeric): length of the coding sequence in Kbp

• oldest_phylostratum (numeric): a number indicating the evolutionary age of the gene (from 1 = gene that appeared before the last common ancestor of all cellular organism to 20 = human or ape-specific gene; thus a smaller number indicates an older gene)

• gc_cds (numeric): GC content of the coding region

• is_ab (logical): whether the gene is known to be involved in abnormal brain morphology when mutated The data frame “gene_te” contains the same information for \~8,000 other genes.

## Problem 1 (20 points)

Using the appropriate regression on the gene_tr data frame answer the following questions:

(a) do longer genes tend to have higher or lower GC content than shorter ones?

(b) is this effect significant?

(c) build a linear model predicting the GC content of a gene from its length and its evolutionary age: what fraction of the GC content variance is explained by this model? What is its mean squared error (MSE)?

(d) based on the previous model, is the CDS length of a gene a predictor of GC content independent of evolutionary age?

(e) use the model built at point (c) to predict the GC content of the genes contained in the data frame gene_te from their length and evolutionary age: what is the MSE of this prediction? Is there evidence of overfitting?

## Problem 2 (10 points)

Still using the gene_tr data frame and the appropriate regression:

(a) which among CDS length, GC content, and evolutionary age is, by itself, a significant predictor of involvement in brain abnormal morphology?

(b) are older genes more or less likely to be involved in brain abnormal morphology? Produce a graph demonstrating your answer to this question

(c) which among CDS length, GC content, and evolutionary age is a significant predictor of involvement in brain abnormal morphology independent of the other two? 1

## Problem 1

### Part (a) - Relationship between gene length and GC content

```{r}
# Load data
load("test_2406214.RData")

# Linear regression to analyze gene length and GC content
lreg <- lm(gc_cds ~ cds_length, data = gene_tr)
summary(lreg)
```

The coefficient for cds_length is negative, indicating that longer genes tend to have lower GC content.

### Part (b) - Significance of the relationship

```{r}
# Check p-value for cds_length
summary(lreg)$coefficients
```

The p-value associated with cds_length is less than 0.05, indicating that the relationship between gene length and GC content is significant.

### Part (c) - Linear model predicting GC content from gene length and evolutionary age

```{r}
# Build linear model including evolutionary age
lreg_1 <- lm(gc_cds ~ cds_length + oldest_phylostratum, data = gene_tr)
summary(lreg_1)
```

The 𝑅\^2 value is approximately 1.39%, indicating a very low fraction of variance explained. The MSE for the full model is 0.007610067.

### Part (d) - Independence of gene length as a predictor of GC content

```{r}
# Compare models and check p-value for cds_length in lreg_1
summary(lreg_1)
```

The p-value for `cds_length` in the full model (`lreg_1`) is 0.0552, suggesting that gene length is not a significant predictor of GC content independent of evolutionary age.

### Part (e) - Predicting GC content for gene_te and evaluating MSE

```{r}
# Predict GC content for gene_te using lreg and lreg_1
pred_gene_te <- predict(lreg, newdata = gene_te)
pred_gene_te1 <- predict(lreg_1, newdata = gene_te)

# Calculate MSE
mse_te <- mean((gene_te$gc_cds - pred_gene_te)^2)
mse_te1 <- mean((gene_te$gc_cds - pred_gene_te1)^2)

mse_te
mse_te1
```

The MSE for the reduced model (lreg) is 0.007555024 and for the full model (lreg_1) is 0.00749084, indicating no evidence of overfitting.

## Problem 2

### Part (a) - Predictors of involvement in abnormal brain morphology

```{r}
# Logistic regression to analyze predictors for abnormal brain morphology
logreg1 <- glm(is_ab ~ cds_length, family = "binomial", data = gene_tr)
logreg2 <- glm(is_ab ~ gc_cds, family = "binomial", data = gene_tr)
logreg3 <- glm(is_ab ~ oldest_phylostratum, family = "binomial", data = gene_tr)

summary(logreg1)
summary(logreg2)
summary(logreg3)
```

Gene length (cds_length) and evolutionary age (oldest_phylostratum) show significant predictors for abnormal brain morphology.

### Part (b) - Relationship between evolutionary age and abnormal brain morphology

```{r}
# Boxplot to visualize relationship between oldest_phylostratum and is_ab
boxplot(oldest_phylostratum ~ is_ab, data = gene_tr)
```

From the plot, older genes are less likely to be involved in abnormal brain morphology.

## Part (c) - Independent predictor of abnormal brain morphology

```{r}
# Logistic regression including all predictors
logreg <- glm(is_ab ~ cds_length + gc_cds + oldest_phylostratum, family = "binomial", data = gene_tr)
summary(logreg)
```

Evolutionary age (oldest_phylostratum) appears to be the most significant predictor, while gene length (cds_length) also shows significant results.
