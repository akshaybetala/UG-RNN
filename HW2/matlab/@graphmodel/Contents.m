function Contents(gm)
% Graphical Model Base Class 
%   represents a basic collection of variables and factors
% 
% METHODS
% Constructors
%   graphmodel( flist ) : construct a graphical model from a collection of factors
%
% Accessors & Factor Manipulation
%   display, disp       : display the basic properties of the model
%   addFactor           : add a factor to the collection
%   getFactor           : access a subset of the factors
%   removeFactor        : remove a factor from the collection
%   withVariable        : get indices of factors that involve a given set of variables
%
% Boolean Properties
%   isbinary(gm)        : are all the variables binary-valued?
%   ispairwise(gm)      : are all the factors pairwise (0-2 variables)?
%
% Inference and Complexity
%   inducedWidth(gm,order) : compute the maximal clique encountered for a given elimination order
%   order(gm,Method)       : find an elimination order using various heuristic methods
%   jtree(gm)              : compute marginals & partition function by junction tree inference
%   joint                  : compute the full joint distribution as a single factor

% TODO: add new functions, increase comments, etc

%   (markovBlanket **???**)
%
% Inference and Complexity
%   inducedWidth(gm,order) : compute the maximal clique encountered for a given elimination order
%   orderMinFill(gm)       : find an elimination order using the min-fill heuristic
%   orderMinWidth(gm)      : find an elimination order using the min-width heuristic
%   orderRandom(gm)        : return a randomly chosen elimination order
% 
%   joint                  : compute the full joint distribution as a single factor
%   MBE*                   : inference using (mini-) bucket elimination
%   DynaDecomp (2,+,+2)    : distribution-approximating elimination procedures
%   
%
% consolidate (?)
% 
  help graphmodel/Contents

