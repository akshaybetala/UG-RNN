clear mex;

%%%--ACCESSOR FUNCTIONS--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mex -I../../cpp/include mexFactorIdx.cpp -DIND2SUBV -output ../ind2subv
mex -I../../cpp/include mexFactorIdx.cpp -DVALUE -output ../value      
mex -I../../cpp/include mexFactorIdx.cpp -DSUBV2IND -output ../subv2ind


%%%--BINARY OPERATORS--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  mex -I../../cpp/include mexFactorBinOp.cpp -DPLUS  -output ../plus
  mex -I../../cpp/include mexFactorBinOp.cpp -DMINUS -output ../minus
  mex -I../../cpp/include mexFactorBinOp.cpp -DTIMES -output ../times
  mex -I../../cpp/include mexFactorBinOp.cpp -DRDIV  -output ../rdivide

%%%--ELIMINATION OPERATORS--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  mex -I../../cpp/include mexFactorElimOp.cpp -DSUM         -output ../sum
  mex -I../../cpp/include mexFactorElimOp.cpp -DLOGSUMEXP   -output ../logsumexp
  mex -I../../cpp/include mexFactorElimOp.cpp -DMAX         -output ../max
  mex -I../../cpp/include mexFactorElimOp.cpp -DMIN         -output ../min
  mex -I../../cpp/include mexFactorElimOp.cpp -DMARGINAL    -output ../marginal
  mex -I../../cpp/include mexFactorElimOp.cpp -DMAXMARGINAL -output ../maxmarginal
  mex -I../../cpp/include mexFactorElimOp.cpp -DMINMARGINAL -output ../minmarginal

  mex -I../../cpp/include mexFactorElimOp.cpp -DCONDITION   -output ../condition

%%%--GROUP OPERATORS--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  mex -I../../cpp/include mexFactorMean.cpp -DMEAN    -output ../mean
  mex -I../../cpp/include mexFactorMean.cpp -DGEOMEAN -output ../geomean

%%%--DISTANCE FUNCTIONS--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  mex -I../../cpp/include mexFactorDist.cpp -DENTROPY -output ../entropy
  %mex mexFactorDist.cpp -DDIST
  %movefile(['mexFactorDist.',mexext],['../distance.',mexext]);				% not quite the same "type" interface!
  %mex mexFactorDist.cpp -DNORM
  %movefile(['mexFactorDist.',mexext],['../norm.',mexext]);					% same problem


%%%--DECOMPOSITION FUNCTIONS--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  mex -I../../cpp/include mexFactorDecomp.cpp -DSUM  -output ../decompSum
  mex -I../../cpp/include mexFactorDecomp.cpp -DPROD -output ../decompProd
