-----------------------------------------------------------------------------
--
-- Module      :  Sobol
-- Copyright   :  (c) Haskellers 2016
-- License     :  BSD-style (see the file LICENSE)
--
-- Maintainer  :  c@cdprf.me
-- Stability   :  development
-- Portability :  -------
--
-- This libraty generates Sobol Sequeneces
--
-- This implementation uses Gray codes, reference:
--
--   * Antonov, I.A. and Saleev, V.M. (1979) "An economic method of
--     computing LPτ-sequences". Zh. Vych. Mat. Mat. Fiz. 19:
--     243–245 (in Russian); U.S.S.R Comput. Maths. Math. Phys. 19:
--     252–256 (in English).
--
-----------------------------------------------------------------------------

module Sobol
  (
  ) where

import Data.Bits

-- SobolState = ( exp, x, idx, poly, polyG, directionalNumbers )
data SobolState = SobolState Int Float Int Int Int [Int]

getPolynomialDegree :: Int -> Int
getPolynomialDegree p = _getPolynomialDegree p 0 0 where
  _getPolynomialDegree 0 _ pos = pos
  _getPolynomialDegree p idx pos = _getPolynomialDegree (shiftR p 1) (idx+1) ( newPos idx pos ( p .&. 1 == 1) ) where
    newPos i p True = i
    newPos i p False = p

getFirstZero :: Int -> Int
getFirstZero n = _getFirstZero n 1 where
  _getFirstZero 0 idx = idx
  _getFirstZero n idx
    | n .&. 1 == 0 = idx
    | otherwise   = _getFirstZero (shiftR n 1) (idx+1)
