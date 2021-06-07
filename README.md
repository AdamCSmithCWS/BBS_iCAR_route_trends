# BBS_iCAR_route_trends

A spatial iCAR model to estimate route-level trends from BBS data.

The model is written in Stan. It is based on the [bbsBayes](https://github.com/BrandonEdwards/bbsBayes), slope model (e.g., [Sauer and Link 2011](https://doi.org/10.1525/auk.2010.09220)), but without the random year-effects, and with no particular stratification. It's effectively an over-dispersed Poisson regression model with random slopes representing the trends at each route, random intercepts representing the mean abundance (mean count) on each route, as well as the among and within observers effects (true observers, see below)

## Observer effects

The observer effects here are random effects for each observer, not the observer-route combinations, used in the bbsBayes models. This is a new thing for the hierarchical Bayesian BBS models. The MCMC algorithm of JAGS and BUGS has a great deal of trouble separately estimating observer effects from route-level intercepts. The HMC algorithm in Stan, plus the spatially explicit estimates of the route-level intercepts, appears to have much more success.

## Example output - Baird's Sparrow 2004-2019

Baird's Sparrow trends in the Great Plains. Dots represent BBS routes on which Baird's Sparrow has been observed between 2004 and 2019. Over this time period, there's an interesting spatial pattern: declines in the Northeast, increases in the Southwest.

![](Baird's_Sparrow_Trends_2004.png)

## Example output - Chipping Sparrow 2004-2019

Chipping Sparrow trends over the last 15 years suggest the species has been increasing in the St-Lawrence River valley, in the southeastern portion of its range, and through some parts of the Great Plains.

![](Chipping_Sparrow_Trends_2004.png)

## Example output - Wood Thrush 2004-2019

Wood Thrush trends over the last 15 years suggest the species has recently increased in the eastern parts of its range, and some portions of the Appalachians, but generally decreased in the regions where it is most abundant.

![](Wood_Thrush_Trends_2004.png)
