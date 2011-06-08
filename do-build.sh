#!/bin/bash

( cd seqloc && \
  runhaskell Setup.hs clean && \
  runhaskell Setup.hs configure --user -p && \
  runhaskell Setup.hs build && \
  runhaskell Setup.hs haddock && \
  runhaskell Setup.hs install ) && \
( cd seqloc-datafiles && \
  runhaskell Setup.hs clean && \
  runhaskell Setup.hs configure --user -p && \
  runhaskell Setup.hs build && \
  runhaskell Setup.hs haddock && \
  runhaskell Setup.hs install )
