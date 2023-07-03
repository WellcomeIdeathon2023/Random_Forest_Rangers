# Wellcome Ideathon Infectious Disease Challenge 1

## useful links
ideathon data:
https://github.com/WellcomeIdeathon2023/infectious_disease_challenges/tree/main/correlates/datasets

our page:
https://github.com/WellcomeIdeathon2023/Random_Forest_Rangers

co-connect docs:
https://co-connect.ac.uk/docs/

## next steps
- set up git repo ✓
- familiarise with git commits/merges ✓
- look at the data files ✓
  - decide how best to merge/analyse the data ←
  - work out how to import into R and merge ←

## data
- SDY180
  - healthy donor vaccinated with:
    - influenza vaccine or
    - pneumococcal vaccine or
    - placebo
  - blood samples taken before and after
  - outcome:
    - gene expression
    - serum antibody results
    - cytokine analysis
  - time period:
    - 2009/2010
- SDY296/301
  - healthy donor vaccinated with:
    - influenza vaccine only
  - blood samples taken before and after
  - outcome:
    - gene expression
    - serum antibody results
  - time period:
    - 2011/2012 (296)
    - 2012/2013 (301)

# thoughts
why do this:
- targetting vulnerable populations of disease
- finding populations that respond poorly to vaccinations
- optimising the cost (money / logistical) of 

## statistical approaches
- traditional:
  - univariate
  - multivariate
  - LASSO etc..
- modern:
  - kmeans
  - pca 
  - louvain / phenograph
  - tsne/umaps
  - neural network stuff
  - random forest

## challenges of this challenge
- missing data
- types of data
  - many tables of connected variables
  - can connect them using graphs
  - but some of these connections may not make sense (i.e. NA columns)
  - need to interpret a column to see whether we can join on that or not
  - not an obvious heuristic to do this
    - factor levels > 1
    - ending in `ACCESSION`
  - many-to-one mappings, and ranking graph
- going forward and using KNN to input new data
- data standardisation

## challenges of this in general
- disparate data systems
- covid:
  - shielding/vulnerable

## questions we need to answer
- 1. working with biostatisticians:
  - importance of statistical methods
  - valid statistics
  - multiple testing
  - validation / training
- 2. new diseases/vaccines:



