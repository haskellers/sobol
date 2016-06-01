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
-- This libraty generates Sobol Sequeneces, these sequences
-- are called quasi-random sequences or low discrepancy
-- sequences. The library will provide methods to create
-- 1 to N dimensional sobol sequences.
--
-- Ideas:
--   - Provide SobolSeq, an instance of the RandomGen. Following the
--     design in System.Random also provide some functions which
--     return a SobolSeq:
--
--       * mkFixedSobolSeq :: Int-> SobolSeq
--         From a integer seed
--         returns a SobolSeq. Later on the code
--         we explain how is the seed process to generate different
--         sobol sequences as any random number generator.
--
--       * mkSobolSeq :: Int -> [Int] -> SobolSeq
--         Given a primitive polynomial over GS(2) represented by an
--         integer and a list of initial numbers returns a SobolSeq
--
--       * mkPolySobolSeq :: RandomGen r => Int -> r -> SobolSeq
--         Given a primitive polynomial over GS(2) represented by an
--         integer and a random generator returns a SobolSeq
--
-- Notes: This implementation uses Gray codes, reference:
--
--   * Antonov, I.A. and Saleev, V.M. (1979) "An economic method of
--     computing LPτ-sequences". Zh. Vych. Mat. Mat. Fiz. 19:
--     243–245 (in Russian); U.S.S.R Comput. Maths. Math. Phys. 19:
--     252–256 (in English).
--
-- Definitions:
--   * Given the primitive polynomial of degree d
--        x^d + a_1 x^d-1 + a_2 x^d-2+.....+ a_d-1 x + 1
--     is over GS(2) if a_k is 0 or 1 for all k in {1..d-1}
--
--
-----------------------------------------------------------------------------

module Sobol
  (
  ) where

import Data.Bits

mkSobolSeq :: Int -> Int
mkSobolSeq x = x

mkSobolSeq :: Int -> [Int] -> Int
mkSobolSeq x y = x

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
