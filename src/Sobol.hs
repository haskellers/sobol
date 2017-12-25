{-# LANGUAGE TemplateHaskell #-}
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
--     The polynomials will be represented by an integer
--     where each bit represent the a_k, so, the index of the
--     last bit set to '1' is the degree of the polynomial.
--     For example:
--      3 => x + 1 ( degree 1 )
--      7 => x^2 + x + 1 ( degree 2 )
--
--   * Initial numbers, they are a sequence of integer [1,3,37]
--
--   * The "seed" or the initial parameters of a Sobol sequence
--     are a pollynomial and a initial numbers which must follow
--     the rules below:
--       - All initial numbers must be odd
--       - b_i < 2^i for all b_i in [b_0,b_1,....,b_n] where n
--         is the degree of the primitive polynomial
--
-----------------------------------------------------------------------------

module Sobol
  (
  SobolSeq,
  random,
  mkSobolSeq
  ) where

import Data.List.Split
import Data.Bits
import Data.FileEmbed
import Paths_sobol

data SobolSeq = SobolSeq SobolState deriving(Show)

mkSobolSeq :: Int -> [Int] -> SobolSeq
mkSobolSeq p n = SobolSeq $
                 SobolState 0 0 0 p pD (initialNumbers pD n) where
  pD = getPolynomialDegree p

random :: SobolSeq -> (Float,SobolSeq)
random (SobolSeq (SobolState exp x idx p pG initN))
  | exp > c = randomStep0 c (SobolState exp x idx p pG initN)
  | exp < c = randomStep1 c (SobolState exp x idx p pG initN)
  | otherwise = randomStep2 c (SobolState exp x idx p pG initN) where
    c = getFirstZero idx

-----------------------------------------------------------------------------
-- Private functions
-----------------------------------------------------------------------------

-- SobolState = ( exp, x, idx, poly, polyG, directionalNumbers )
data SobolState = SobolState Int Int Int Int Int [Int]
  deriving (Show,Eq)

-- Polynomias
readFile getDataFileName "polynomials.txt"
getPolynomials = parse $ (map (splitOn ",") . (splitOn "\n") $ ($(embedStringFile $ getDataFileName "polynomials.txt") :: String) ) where
  parse = (map $ (map (read :: [Char] -> Int) ))


-- These three methods are which generates the random sequences
randomStep0 :: Int -> SobolState -> (Float,SobolSeq)
randomStep0 c (SobolState exp x idx p pG initN) =
  (r_x, SobolSeq $ SobolState n_exp n_x (idx+1) p pG n_initN) where
    n_initN = computeDirectionalNumbers p pG c initN
    n_exp = exp
    r_x = (fromIntegral n_x ::Float)  / (fromIntegral $ shiftL (1::Int) n_exp :: Float)
    n_x = xor xN xN1 where
      xN = x
      xN1 = n_initN !! (c - 1) * (shiftL (1::Int) (exp - c))

randomStep1 :: Int -> SobolState -> (Float,SobolSeq)
randomStep1 c (SobolState exp x idx p pG initN) =
  (r_x, SobolSeq $ SobolState n_exp n_x (idx+1) p pG n_initN) where
    n_initN = computeDirectionalNumbers p pG c initN
    n_exp = c
    r_x = (fromIntegral n_x :: Float) / (fromIntegral $ shiftL (1::Int) n_exp :: Float)
    n_x = xor xN xN1 where
      xN = x * shiftL (1::Int) (c - exp)
      xN1 = n_initN !! (c - 1)

randomStep2 :: Int -> SobolState -> (Float,SobolSeq)
randomStep2 c (SobolState exp x idx p pG initN) =
  (r_x, SobolSeq $ SobolState n_exp n_x (idx+1) p pG n_initN) where
    n_initN = computeDirectionalNumbers p pG c initN
    n_exp = exp
    r_x = (fromIntegral n_x :: Float) / (fromIntegral $ shiftL (1::Int) n_exp :: Float)
    n_x = xor xN xN1 where
      xN = x
      xN1 = n_initN !! (c - 1)

-- Check and return a inital numbers for a polynomial of
-- polyDegree.
-- The properties that must satisfies are specified on
-- the head of the file.
initialNumbers :: Int -> [Int] -> [Int]
initialNumbers 1 _ = [1]
initialNumbers polyDegree initN
    | any id $ map (\(x,y)->even x || x >= y)
             $ zip initN [ 2^i | i<-[1..(length initN)]]
               = fail "Provide initial numbers that satisfies the properties"
    | length initN >= polyDegree = take polyDegree initN
    | otherwise = fail "Please, provice enough initial numbers"

-- Generate the list of size idx of the directional numbers
-- for the polynomial p and the initial numbers init
-- take pG initN
computeDirectionalNumbers :: Int -> Int -> Int -> [Int] -> [Int]
computeDirectionalNumbers 0 _ idx _ = take idx [1,1..]
computeDirectionalNumbers p pG idx initN
  | (idx - 1) < q = initN
  | otherwise = computeDirectionalNumbers p pG idx $ initN ++  [computeDirectionalNumber p pG newInit] where
    q = length initN
    newInit = initN ++ [initN !! (q - pG)]
computeDirectionalNumber :: Int -> Int -> [Int] -> Int
computeDirectionalNumber p pG initN = foldl xor (last initN) $
                                                zipWith (*) dirNums $
                                                  zipWith (*) powdos polyseq where
  dirNums = drop (length initN - pG) initN
  --powdos = [2^i | i <- [ length dirNums - j | j <- [0..(length dirNums - 1)]]]
  powdos = [2^i | i <- [ pG - j | j <- [0..(pG - 1)]]]
  polyseq = [ shiftR p i .&. 1 | i <- [0..pG-1]]


-- Given a primitive polynomial over GS(2) represented by an
-- integer returns the polynomial degree
-- Example:
--      getPolynomialDegree 3 = 1 (x + 1)
--      getPolynomialDegree 7 = 2 (x^2 + x + 1)
getPolynomialDegree :: Int -> Int
getPolynomialDegree p = _getPolynomialDegree p 0 0 where
  _getPolynomialDegree 0 _ pos = pos
  _getPolynomialDegree p idx pos =
    _getPolynomialDegree (shiftR p 1) (idx+1) ( newPos idx pos ( p .&. 1 == 1) ) where
      newPos i p True = i
      newPos i p False = p

-- Returns the index of the first zero bit from the rigth
-- Examples:
--    getFirstZero 1 = 2 (001)
--    getFirstZero 3 = 3 (011)
getFirstZero :: Int -> Int
getFirstZero n = _getFirstZero n 1 where
  _getFirstZero 0 idx = idx
  _getFirstZero n idx
    | n .&. 1 == 0 = idx
    | otherwise   = _getFirstZero (shiftR n 1) (idx+1)
