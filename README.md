# Code for performing \betaARMA 

## AUTHOR:

- Based on the original code by Fabio M. Bayer (bayer@ufsm.br) on 2015-10-15
- Modified and improved by Everton da Costa (everton.ecosta@ufpe.br) on 2025-04-02
- Modified by Beatriz Ariadna da Silva Ciriaco (beatriz.ciriaco@ufpe.br) on 2025-14-07

## DESCRIPTION:

The code below was originally developed by Fabio M. Bayer (bayer@ufsm.br), 2015-10-15. It includes modifications and improvements made by Everton da Costa (everton.ecosta@ufpe.br) and parcially modified later on by Beatriz Ariadna da Silva Ciriaco only for the BARMA model. 

This code implements beta-ARMA conditional maximum likelihood estimation using in the predictor scale and the original scale for BARMA. It was developed by Fabio M. Bayer, modified by Everton da Costa, and complemented by Beatriz Ciriaco. The modifications include (but are not restricted to): 

- Rewrote the code partially
- Reviewed the code against the original article:
   Rocha, A.V.; Cribari-Neto, F. (2009) Beta autoregressive moving average 
   models, TEST, 18(3), 529-545 (Erratum: 26, 2017, 451-459).
- Optimized the code to avoid computing the same quantity several times 
