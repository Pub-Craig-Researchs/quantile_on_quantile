function [w, sig, c, c1, k] = bds (series, maxdim, distance, flag, maxram)

%BDS Brock, Dechert & Scheinkman test for independence based on the correlation dimension
%
% [W, SIG, C, C1, K] = BDS (SERIES, MAXDIM, DISTANCE, METHOD, MAXRAM)
%
% uses       - time-series vector SERIES (1),
%            - dimensional distance
%              * either defined as fraction DISTANCE of the standard deviation of SERIES
%                if FLAG = 0,
%              * or defined such that the one dimensional correlation integral of SERIES
%                is equal to DISTANCE if FLAG = 1 (2),
%            - not more than MAXRAM megabytes of memory for the computation (3),
%
% to compute - BDS statistics W for each dimension between 2 and MAXDIM (4),
%            - significance levels SIG at which the null hypotheses of no dependence are
%              rejected ASYMPTOTICALLY (use companion function BDSSIG.M for finite
%              samples) against (almost) any type of linear and non-linear dependence (5), 
%            - correlation integral estimates C for each dimension M between 2 and MAXDIM,
%            - first-order correlation integral estimates C1 computed over the last N-M+1
%              observations, and
%            - parameter estimate K (6).
%
% (1) SERIES is normally a vector of residuals obtained from a regression, but it can also
%     be any other stationary time series.
% (2) The default settings are DISTANCE = 1.5 and FLAG = 0. The BDS statistic appears to
%     be most efficient estimated if the measure of dimensional distance EPSILON is chosen
%     such that the first-order correlation integral estimate (C1) lies around 0.7 (see
%     Kanzler, 1998, forthcoming). For settings DISTANCE = 0.7 and FLAG = 1, the
%     programme will chose EPSILON accordingly. Unfortunately, the cost of finding optimal
%     EPSILON is quite high in terms of CPU time and required memory. For a near-normal
%     distribution, the default settings achieve the same without any extra computational
%     burden.
% (3) The default setting is MAXRAM = 150, which is recommended for a system with 192MB
%     physical RAM installed. The programme is highly optimised as to maximise speed given
%     available memory, so it is very important to specify MAXRAM correctly as the amount
%     of physical memory available AFTER starting MATLAB, loading any data and running
%     other applications concurrently. The smaller the amount of RAM available to the
%     programme (in relation to the length of SERIES), the slower the algorithm chosen
%     from six alternatives. However, if MAXRAM is chosen too large, MATLAB will make use
%     of virtual (hard-disk) memory, and this will slow down computation considerably.
% (4) The default setting for MAXDIM is 2. For MAXDIM = 1, W and SIG are empty.
% (5) A vector of NaN is returned if the MATLAB Statistics Toolbox is not installed.
% (6) The BDS statistic W(M) is a function of C(M), C1(1), C1(M) and K, and these
%     estimates are normally of no further interest.
%
% See Kanzler (1998) for some explanation of the main parts of the algorithm (other
% explanations are commented into the below code), for a detailed investigation of the
% finite-sample properties of the BDS statistic, for tables of small-sample quantiles and
% for a comparison with software by Dechert (1988) and LeBaron (1988, 1990, 1997a, 1997b).
% These and other important references can be found at the end of the script.
%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% * All rights reserved. This script may be redistributed if it is left unaltered in    *
% * its entirety (619 lines, 31422 bytes) and if nothing is charged for redistribution. *
% * Usage of the programme in applications and alterations of the code should be        *
% * referenced properly. See http://users.ox.ac.uk/~econlrk for updated versions.       *
% * The author appreciates suggestions for improvement or other feedback.               *
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%
%                    Copyright (c) 14 April 1998  by Ludwig Kanzler
%                    Department of Economics,  University of Oxford
%                    Postal: Christ Church, Oxford OX1 1DP, England
%                    E-mail:  ludwig.kanzler@economics.oxford.ac.uk
%                    $ Revision: 2.41 $ $ Date: 15 September 1998 $

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % Executable part of main function BDS.M starts here  % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%%%%%%%%%%%%%%%%%%%%%% Check and transformation of input arguments %%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
   maxram = 150;
elseif maxram > 500
   disp('Are you sure you have so much memory available?')
   error('If so, you need to edit the code, otherwise try again with a lower value.')
end

if nargin < 4
   flag = 0;
elseif ~any(flag == [0 1])
   error('Unknown method for determining dimensional distance; try again with 0 or 1.')
end

if nargin < 3
   distance = 1.5;
elseif distance < 0
   error('The dimensional distance parameter must be positive.')
elseif flag == 1 & distance > 1
   error('The correlation integral cannot exceed 1.')
end

if nargin < 2
   maxdim = 2;
elseif maxdim < 1
   error('The dimension needs to be at least 1.');
end

if nargin < 1
   error('Cannot compute the BDS statistic on nothing.')
end

[rows,cols] = size(series);
if rows > 1 & cols == 1
   n = rows;
   series = series';
elseif cols > 1 & rows == 1
   n = cols;
elseif cols > 1 & rows > 1
   n = cols*rows;
   series = series(:)'; % transformation into a row vector
   disp(sprintf('\aTransformed matrix input into a single column.'))
else
   error('Cannot compute the BDS statistic on a scalar!')
end

%%%%%%%%%%%% Determination of and preparations for fastest method given MAXRAM %%%%%%%%%%%

fastbuild     = 0.000016 * (1:52) .* pow2(1:52); % memory requirements
slowbuild     = 0.000045           * pow2(1:52); % for the various
holdinfo      = 0.000005           * pow2(1:52); % algorithms in 
wordtable     = 0.000008 *   n^2  ./     (1:52); % megabytes for 
bitandop      = 0.000024 *   n^2  ./     (1:52); % given N

[ram1, bits1] = min(fastbuild + holdinfo + wordtable + bitandop); % number of bits for
[ram2, bits2] = min(fastbuild + holdinfo + wordtable);            % which each of six
[ram3, bits3] = min(slowbuild + holdinfo + wordtable + bitandop); % methods uses minimum
[ram4, bits4] = min(slowbuild + holdinfo + wordtable);            % memory; this memory
[ram5, bits5] = min(                       wordtable + bitandop); % is given by
[ram6, bits6] = min(                       wordtable);            % ram1, ram2,..., ram6

if ram1 < maxram | ram2 < maxram
   if ram1 < maxram
      method = 1;
      bits = bits1; ram = ram1;
   else
      method = 2;                                     % maximum number of rows to put
      bits = bits2; ram = ram2;                       % through BITAND and bit-counting
      stepping = floor((maxram-ram)*bits/n/0.000024); % algorithm without exceeding MAXRAM
   end

   % Vector BITINFO lists the number of bits set for each integer between 0 and 2^bits
   % (corresponding to the indices of the vector shifted by 1). See Kanzler (1998) for an
   % explanation.
   bitinfo = uint8(sum(rem(floor((0:pow2(bits)-1)'*pow2(1-bits:0)),2),2));

elseif ram3 < maxram | ram4 < maxram
   if ram3 < maxram
      method = 3;
      bits = bits3; ram = ram3;
   else
      method = 4;
      bits = bits4; ram = ram4;
      stepping = floor((maxram - ram) * bits / n / 0.000024);
   end

   bitinfo(1:pow2(bits), :) = uint8(0);         % the same as above, but created through
   for bit = 1 : bits                           % a loop, which consumes less memory
      bitinfo(1:pow2(bits)) = sum([bitinfo, ...
         kron(ones(pow2(bits-bit),1), [zeros(pow2(bit-1),1); ones(pow2(bit-1),1)])],2);
   end

elseif ram5 < maxram | ram6 < maxram
   if ram5 < maxram
      method = 5;
      bits = bits5; ram = ram5;
   else
      method = 6;
      bits = bits6; ram = ram6;
      stepping = floor((maxram - ram) * bits / n / 0.000024);
   end

else
   disp('Insufficient amount of memory. Allocate more memory to the system')
   disp('or reduce the number of observations, then try again.')
   error(' ')
end

%%%%%%%%%%%%%%%%%%%%% Determination of dimensional distance EPSILON %%%%%%%%%%%%%%%%%%%%%%

% The empirical investigation by Kanzler (1998) shows that choosing EPSILON such that the
% first-order correlation integral is around 0.7 yields the most efficient estimation of
% low-dimensional BDS statistics. Hence the objective here is to choose EPSILON such that,
% say, 70% of all observations lie within distance EPSILON to each other. If desired, the
% programme first determines EPSILON as to fulfil this or a similar requirement.
%
% The conceptually simplest way of setting up the calculation of distance among all
% observations is to define a two-dimensional table D (for "distance") of length and width
% N and assign to each co-ordinate (x,y) the result of the problem ABS(x-y).
%
% In principle, the entire table could thus be created with the following one-line
% statement:
%             D = ABS( SERIES(ONES(1000,1),:)' - SERIES(ONES(1000,1),:) )
%
% Since the lower triangle of the table only replicates the upper triangle and since the
% diagonal values represent own values (ones) which are not desired to be included in the
% calculation, only the upper triangle receives further attention.
%
% Unfortunately, sewing all the row vectors of the upper triangle together to form one
% single (row) vector makes indexing very messy. To aid understanding of the vector-space
% indexing used here (as well as in the optional sub-function further below), one may wish
% to refer to the following exemplary matrix table (N=7):
%
%                                       Using this example, it is easy to verify
%          * * * *c o l u m n* * * *    that column vector I is defined by the
%     I    1   2   3   4   5   6   7    following indices in vector space:
%                                          I+(0 : I-2)*N - CUMSUM(1 : I-1)
%  *  1    *   1   2   3   4   5   6
%  *  2    .   *   7   8   9  10  11    More generally, column vector I starting only
%  r  3    .   .   *  12  13  14  15    in row J is:
%  o  4    .   .   .   *  16  17  18       I+(J-1 : I-2)*N - SUM(1:J-1)-CUMSUM(J : I-1)
%  w  5    .   .   .   .   *  19  20
%  *  6    .   .   .   .   .   *  21    Row vector I is given by indices:
%  *  7    .   .   .   .   .   .   *       1+(I-1)*(N-1)-SUM(1:I-2) : I*(N-1)-SUM(1:I-1)
%
% (A formal derivation of the above formulae is beyond the scope of this script.)
%
% To calculate a percentile of the distribution of distance values, the row vector is
% sorted (unfortunately, this requires a lot of time and RAM in MATLAB).

if ~flag
   demeaned = series-sum(series)/n;                        % fastest algorithm for
   epsilon  = distance * sqrt(demeaned*demeaned'/(n-1));   % computing the standard
   clear demeaned % to save memory                         % deviation of SERIES
   
elseif 0.000008 * 3 * sum(1:n-1) < maxram % check memory requirements for DIST and sorting
   dist(1:sum(1:n-1)) = 0;
   for i = 1 : n-1
      dist(1+(i-1)*(n-1)-sum(0:i-2):i*(n-1)-sum(1:i-1)) = abs(series(i+1:n)-series(i));
   end
   sorted  = sort(dist);
   epsilon = sorted(round(distance*sum(1:n-1))); % DISTANCEth percentile of SORTED series
   clear dist sorted
else
   error('Insufficient RAM to compute EPSILON; allocate more memory or use METHOD = 1.')
end

%%%%%%%%%%%% Computation and storage of one-dimensional distance information %%%%%%%%%%%%

% Similarly to the above, a two-dimensional table C (for "close") of length and width N
% can be defined by assigning to each co-ordinate (x,y) the result of the problem
% ABS(x-y) <= EPSILON; (x,y) assumes the value 1 if the statement is true and 0 otherwise.
%
% Formally, for given EPSILON:
%                                 C(x,y) = 1 if ABS(x-y) <= EPSILON
%                                        = 0 otherwise
%
% Once again, the resulting information needs to be stored in the most efficient way.
% In this implementation, this is done by chopping each row of the table into "words" of
% several bits, the precise number of bits per word being determined by the above
% algorithms. One "word" is thus represented by one integer. This slashes the size of the
% table by the number of bits. See Kanzler (1998) for more details.
%
% The below routine stores all rows of the upper triangle of the conceptual table
% (described in Kanzler, 1998) left-aligned and assigns zeros to all other elements.
%
% As will also be explained further below, the computation of parameter K requires the sum
% of each FULL row, i.e. each row including the elements in the lower triangle and on the
% diagonal. The "missing" bits correspond to the sums over each column in the upper
% triangle, and these sums are also computed and stored in the below loop. And to make
% matters simple, diagonal values are allocated to the column sums by initialising them
% with value 1. See also Kanzler (1998).

colsum(1:n)              = 1;
rowsum(1:n)              = 0;
nwords                   = ceil((n-1)/bits);
wrdmtrx(1:n-1,1:nwords)  = 0;                 % initialisation of bit-word table

for row = 1 : n-1
   bitvec                = abs(series(1+row:n) - series(row)) <= epsilon;
   rowsum(row)           = sum(bitvec);
   colsum(1+row:n)       = colsum(1+row:n) + bitvec;
   nwords                = ceil((n-row)/bits);
   wrdmtrx(row,1:nwords) = (reshape([bitvec,zeros(1,nwords*bits-n+row)],...%transformation
                                        bits, nwords)' *pow2(0:bits-1)')'; %into bit-words
end
clear series bitvec

%%%%%%%%%%%%%%%%%% Computation of one-dimensional correlation estimates %%%%%%%%%%%%%%%%%%

% C1(1), the fraction (or estimated probability) of pairs in SERIES being "close" in the
% first dimension is just the average over ALL unique elements. C1(1) is hence the most
% efficient estimator of C(1), and the resulting estimate is used in the computation of
% SIGMA(M) further below.
%
% However, for the difference term C(M) - C(1)^M of the BDS statistic (see further below)
% to follow SIGMA asymptotically, both C(M) and C(1) need to be estimated over the same
% length vector, and so MAXDIM different C1's need to be estimated here:
%
%                               N     N
%    C1(M) = 2/(N-M+1)/(N-M) * SUM   SUM  B(S,T)
%                              S=M  T=S+1
%
% Each C1(M) is easily computed from the sum of all bits set in rows M to N-1 divided by
% the appropriate total number of bits.

bitsum(maxdim:-1:1) = cumsum([sum(rowsum(maxdim:n-1)), rowsum(maxdim-1:-1:1)]);
c1    (maxdim:-1:1) = bitsum(maxdim:-1:1) ./ cumsum([sum(1:n-maxdim), n-maxdim+1 : n-1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Computation of parameter K %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A parameter needed to estimate SIGMA(M) is K, which is defined as:
%
%                       N     N     N
% K = 6/N/(N-1)/(N-2)* SUM   SUM   SUM {C(T,S)*C(S,R) + C(T,R)*C(R,S) + C(S,T)*C(T,R)} / 3
%                      T=1  T=T+1 R=S+1
%
% As is readily apparent, a literary computation of the above would be very processing
% intensive, e.g.:
%                  HT(1 : N) = 0;
%                  FOR T = 1 : N
%                     HS(1 : N-T) = 0;
%                     FOR S = T+1 : N
%                        HR(1 : N-S) = 0;
%                        FOR R = S + 1 : N
%                           HR(R-S) = (C(T,S)*C(S,R) + C(T,R)*C(R,S) + C(S,T)*C(T,R)) / 3;
%                        END
%                        HS(S-T) = SUM(HR(1 : N-S));
%                     END
%                     HT(T)= SUM(HS(1 : N-T));
%                  END
%                  K = SUM(HT) * 6 / N/(N-1)/(N-2);
%
% To understand what k actually estimates, and how this estimation can be made
% computationally more efficient, see Kanzler (1998).
%
% The above FOR loop computes the sum over each row and over each column including the
% diagonal in the upper triangle. To compute K from this, the sum of the squares of the
% row and column sums needs to be adjusted as reasoned above, whereby the sum of all
% elements in table C is given by twice the sum of all vector elements plus the diagonal
% values.

fullsum = rowsum + colsum;
k       = (fullsum*fullsum' + 2*n - 3*(2*bitsum(1)+n)) / n/(n-1)/(n-2);
clear rowsum colsum fullsum bitsum

%%%%%%%%%% Computation of correlation estimates and SIGMA for higher dimensions %%%%%%%%%%

% C(M), the M-dimension correlation
% estimate, is defined as:                            N     N    M-1
%                           C(M) = 2/(N-M+1)/(N-M) * SUM   SUM   PROD B(S-J, T-J)
%                                                    S=M  T=S+1  J=0
%
% To see how C can be computed for M > 1, see Kanzler (1998).
%
% In practice, the required BITAND-operation can be performed on the entire table at once
% by replacing the entire table between rows M and N-1 with the result of the BITAND-
% operation between the table formed by rows M to N-1 and M-1 to N-2. But this works only
% if sufficient memory is available (methods 1, 3 and 5). Otherwise, the BITAND-operation
% has to be performed by looping BACKWARDS through the table, taking as many rows as
% possible at once (methods 2, 4 and 6).
%
% The number of bits set in rows M to N-1 (inclusive) is counted either in one go or
% through the above loop by looking up the number of bits set for each integer (LeBaron,
% 1997, uses a similar method), or, if memory was insufficient to create the required
% BITINFO array, by column-wise brute-force counting.
%
% MATLAB uses logarithms to compute powers, and this can result in minute deficiencies in
% accuracy. To avoid this, integer powers are computed by separate functions in this
% script (see further below). Otherwise, SIGMA would be calculated as follows:
%    sigma(m-1) = 2*sqrt(k^m + 2*k.^(m-(1:m-1))*(c1(1).^(2*(1:m-1)))'...
%                        + (m-1)^2*c1(1)^(2*m) - m^2*k*c1(1)^(2*m-2));

for m = 2 : maxdim
   bitcount = 0;

   if sum(method == [1 3])
      wrdmtrx(m:n-1,:) = bitand(wrdmtrx(m:n-1,:),wrdmtrx(m-1:n-2,:)); % BITAND and bit
      bitcount = sum(sum(bitinfo(wrdmtrx(m:n-1,:)+1)));               % count all at once

   elseif sum(method == [2 4])
      for row = n-stepping : -stepping : m+1                                   % BITAND
         wrdmtrx(row:row+stepping-1,:) = bitand(wrdmtrx(row:...                % and bit
                          row+stepping-1,:), wrdmtrx(row-1:row+stepping-2,:)); % count in
         bitcount=bitcount+sum(sum(bitinfo(wrdmtrx(row:row+stepping-1,:)+1))); % backward
      end                                                                      % loops
      wrdmtrx(m:row-1,:)  = bitand(wrdmtrx(m:row-1,:), wrdmtrx(m-1:row-2,:));  % through
      bitcount = bitcount + sum(sum(bitinfo(wrdmtrx(m:row-1,:)+1)));           % the table

   elseif method == 5
      wrdmtrx(m:n-1,:) = bitand(wrdmtrx(m:n-1,:), wrdmtrx(m-1:n-2,:)); % BITAND at once...
      for col = 1 : ceil((n-1)/bits)                                   % bit count
         bitcount = bitcount + sum(sum(rem(floor(wrdmtrx(m:...         % by brute force
                       n-1-(col-1)*bits, col) * pow2(1-bits:0)), 2))); % in loops
      end

   else
      for row = n-stepping : -stepping : m+1
         wrdmtrx(row:row+stepping-1,:) = bitand(wrdmtrx(row:...               % BITAND
                   row+stepping-1,:), wrdmtrx(row-1:row+stepping-2,:));       % opera-
      end                                                                     % tions
      wrdmtrx(m:row-1,:) = bitand(wrdmtrx(m:row-1,:), wrdmtrx(m-1:row-2,:));  % and brute-
      for col = 1 : ceil((n-1)/bits)                                          % force bit
         bitcount = bitcount + sum(sum(rem(floor(wrdmtrx(m:...                % counting
            n-1-(col-1)*bits, col) * pow2(1-bits:0)),2)));                    % in loops
      end
   end

   c(m-1)     = bitcount / sum(1:n-m);                                       % indexing of
   sigma(m-1) = 2*sqrt(prod(ones(1,m)*k) + 2*ivp(k,m-(1:m-1),m-1)...         % C and SIGMA
                *(ivp(c1(1),2*(1:m-1),m-1))' + (m-1)*(m-1)...                % runs from 1
                *prod(ones(1,2*m)*c1(1)) - m*m*k*prod(ones(1,2*m-2)*c1(1))); % to MAXDIM-1
end
clear wrdmtrx

%%%%%%%%%%%%%%% Computation of the BDS statistic and level of significance %%%%%%%%%%%%%%%

% Under the null hypothesis of independence, it is obvious that the time-series process
% has the property C(1)^M = C(M). In finite samples, C(1) and C(M) are consistently
% estimated by C1(M) and C(M) as above. Also, Brock et al. (1996) show that the standard
% deviation of the difference C(M) - C1(M)^M can be consistently estimated by SIGMA(M)
% divided by SQRT(N-M+1), where:
%
%                           M-1
% SIGMA(M)^2 = 4* [K^M + 2* SUM {K^(M-J)* C^(2*J)} + (M-1)^2* C^(2*M) - M^2* K* C^(2*M-2)]
%                           J=1
%
% and C = C1(1) and K as above.
%
% For given N and EPSILON, the BDS Statistic                        C(M) - C1(M)^M
% is defined as the ratio of the two terms:    W(M) = SQRT(N-M+1) * --------------
%                                                                      SIGMA(M)
%
% Since it follows asymptotically the normal distribution with mean 0 and variance 1,
% hypothesis testing is straightforward. If available, this is done here using function
% NORMCDF of the MATLAB Statistics Toolbox.
%
% Integer powers are again calculated by a sub-routine which is more accurate than the
% MATLAB built-in power function; without using the sub-routine, the line for calculating
% W would be: w = sqrt(n-(2:maxdim)+1) .* (c - c1(2:maxdim).^(2:maxdim)) ./ sigma;

if maxdim > 1
   w = sqrt(n-(2:maxdim)+1) .* (c - idvp(c1(2:maxdim), 2:maxdim, maxdim-1)) ./ sigma;

   if exist('normcdf.m','file') & nargout > 1
      sig = min(normcdf(w,0,1), 1-normcdf(w,0,1)) * 2;
   elseif nargout > 1
      sig(1:maxdim-1) = NaN;
   end

else
   w   = [];
   sig = [];
   c   = [];
end

%%%%%%%%%%%%%%%%%%%%%%% Sub-functions for computing integer powers %%%%%%%%%%%%%%%%%%%%%%%

function ipow = ivp (base, intpowvec, veclen)
   ipow(1 : veclen) = 0;
   for j = 1 : veclen
      ipow(j) = prod(ones(1, intpowvec(j)) * base);
   end

function ipow = idvp (basevec, intpowvec, veclen)
   ipow(1 : veclen) = 0;
   for j = 1 : veclen
      ipow(j) = prod(ones(1, intpowvec(j)) * basevec(j));
   end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % Executable part of main function BDS.M ends here  % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                                       %
%   The following sub-function is not actually used by the main function and only       %
%   included for the benefit of those who would like to implement the BDS test in a     %
%   language which is either incapable of or inefficient in handling bit-wise AND-      %
%   operations, or those who would like to cross-check the above computation. Deleting  %
%   the sub-function from the script will NOT result in any increase in performance.    %
%                                                                                       %
%   To use the function, save the remainder of this code in a file named BDSNOBIT.M.    %
%                                                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function w = bdsnobit (series, maxdim, eps)
%BDSNOBIT BDS test for independence IMPLEMENTED WITHOUT USING BIT-WISE FUNCTIONS
%
% Only such comments which relate exclusively to this implementation of the test and
% which cannot be found in the main function are included below.
%
%                    Copyright (c) 14 April 1998  by Ludwig Kanzler
%                    Department of Economics,  University of Oxford
%                    Postal: Christ Church, Oxford OX1 1DP, England
%                    E-mail:  ludwig.kanzler@economics.oxford.ac.uk
%                    $ Revision: 1.3 $      $ Date: 30 April 1998 $

%%%%%%%%%%%%%%%%%%%%%% Check and transformation of input arguments %%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
   eps = 1;
   if nargin == 1
      maxdim = 2;
   elseif maxdim < 2
      error('MAXDIM needs to be at least 2!');
   end
end

epsilon = std(series)*eps;
series  = series(:)';      % PAIRS is the total number of unique pairs which can be
n       = length(series);  % formed from all observations (note that while this is 
pairs   = sum(1:n-1);      % just (N-1)*N/2, MATLAB computes SUM(1:N-1) twice as fast!)

%%%%%%%%%%%% Computation and storage of one-dimensional distance information %%%%%%%%%%%%

% Recall that in the implementation of the main function above, table C is stored in bit-
% representation. When this is not possible or desirable, the second best method is to use
% one continuous vector of unassigned 8-bit integers (called UINT8). This, however,
% requires version 5.1 or higher, and a similar option may not be available in other high-
% level languages. Implementation does not depend on the ability to use unassigned low-bit
% integers and would work equally with double-precision integers, but the memory
% requirements would, of course, be higher. Using UINT8's is still a rather inefficient
% way of storing zeros and ones, which in principle require only a single bit each. On the
% PC, MATLAB actually requires "only" around 5 bytes for each UNIT8.

b(1:pairs) = uint8(0);
for i = 1 : n-1
   b(1+(i-1)*(n-1)-sum(0:i-2):i*(n-1)-sum(1:i-1)) = abs(series(i+1:n)-series(i))<=epsilon;
end
clear series

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Computation of parameter K %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sums(1 : n) = 0;
for i = 1 : n
   sums(i) =   sum(b(i+(0 : i-2)*n - cumsum(1 : i-1)))...             % sum over column I
             + 1 ...                                                  % diagonal element
             + sum(b(1+(i-1)*(n-1)-sum(1:i-2) : i*(n-1)-sum(1:i-1))); % sum over row I
end
k = (sum(sums.^2) + 2*n - 3*(2*sum(b)+n)) / n/(n-1)/(n-2);

%%%%%%%%%%%%%%%%%% Computation of one-dimensional correlation estimates %%%%%%%%%%%%%%%%%%

bitsum(1:maxdim) = sum(b(1+(maxdim-1)*(n-1)-sum(0:maxdim-2) : pairs));
for m = maxdim-1 : -1 : 1
   bitsum(m) = bitsum(m+1) + sum(b(1+(m-1)*(n-1)-sum(0:m-2):m*(n-1)-sum(1:m-1)));
end
c1(maxdim:-1:1) = bitsum(maxdim:-1:1) ./ cumsum([sum(1:n-maxdim), n-maxdim+1 : n-1]);

%%%%%%%%%% Computation of correlation estimates and SIGMA for higher dimensions %%%%%%%%%%

for m = 2 : maxdim

   % Indexing in vector space once again follows the rules set out above. Multiplication
   % is done by moving up column by column into  north-west direction, so counter I runs
   % backwards in the below WHILE loop until the Mth column (from the left) is reached:
   i = n;
   while i - m

      % Multiplication is not defined on UINT8 variables and translating the columns
      % twice, once from UINT8 to DOUBLE integer and then back to UINT8, would be
      % inefficient, so it is better to sum entries (this operation - undocumented by
      % MATLAB - is defined, and even faster than the documented FIND function!) and
      % compare them against the value 2:
      b(i + (m-1 : i-2)*n - sum(1:m-1) - cumsum(m : i-1)) = ...
         sum([  b(i   + (m-1 : i-2)*n - sum(1:m-1) - cumsum(m   : i-1)); ...
                b(i-1 + (m-2 : i-3)*n - sum(1:m-2) - cumsum(m-1 : i-2))  ]) == 2;
      
      % The sum over each column is computed immediately after that column has been
      % updated. To store the column sums, the vector SUMS already used above for the row
      % sums is recycled (this is more memory-efficient than clearing the above SUMS
      % vector and defining a new vector of the column sums, because in the latter case,
      % MATLAB's memory space will end up being fragmented by variables K and C added to
      % the memory in the meantime!):
      sums(i) = sum(b(i + (m-1 : i-2)*n - sum(1:m-1) - cumsum(m : i-1)));
   i = i - 1;
   end

   c(m-1)     = sum(sums(m+1:n)) / sum(1:n-m);
   sigma(m-1) = 2*sqrt(k^m + 2*k.^(m-(1:m-1))*(c1(1).^(2*(1:m-1)))'... % could use above
                         + (m-1)^2*c1(1)^(2*m) - m^2*k*c1(1)^(2*m-2)); % inter-power sub-
end                                                                    % functions instead

%%%%%%%%%%%%%%% Computation of the BDS statistic and level of significance %%%%%%%%%%%%%%%

w = sqrt(n-(2:maxdim)+1) .* (c-c1(2:maxdim).^(2:maxdim)) ./ sigma; % or use sub-functions

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % Sub-function BDSNOBIT.M ends here % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


% REFERENCES:
%
% Brock, William, Davis Dechert & Jos� Scheinkman (1987), "A Test for Independence
%    Based on the Correlation Dimension", University of Wisconsin-Madison, Social
%    Science Research Working Paper, no. 8762
%
% Brock, William, Davis Dechert, Jos� Scheinkman & Blake LeBaron (1996), "A test for
%    independence based on the correlation dimension", Econometric Reviews, vol. 15,
%    no. 3 (August), pp. 197-235, revised version of Brock et al. (1987)
%
% Dechert, Davis (1988), "BDS STATS: A Program to Calculate the Statistics of the
%    Grassberger-Procaccia Correlation Dimension Based on the Paper "A Test for
%    Independence" by W. A. Brock, W. D. Dechert and J. A. Scheinkman", version 8.21
%    (latest), MS-DOS software available on gopher.econ.wisc.edu
%
% Kanzler, Ludwig (1998), "Very Fast and Correctly Sized Estimation of the BDS Statistic",
%    Oxford University, Department of Economics, working paper, available on
%    http://users.ox.ac.uk/~econlrk
%
% LeBaron, Blake (1988, 1990, 1997a), "BDSTEST.C", version June 1997 (latest), C source
%    code available on gopher.econ.wisc.edu
%
% LeBaron, Blake (1997b), "A Fast Algorithm for the BDS Statistic", Studies in
%    Nonlinear Dynamics and Econometrics, vol. 2, pp. 53-59


% ACKNOWLEDGEMENT:
%
% I am grateful to Blake LeBaron for giving me the exclusive opportunity to beta-test his
% C programme in its compiled version for MATLAB 5 and thus enabling me to compare the two
% programmes directly. I have benefited from the many associated discussions.


% End of file.
