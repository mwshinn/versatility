% find_nodal_versatility will calculate the versatility of each node
% in a network.  `g` is an NxN adjacency matrix, which can either be
% weighted or unweighted, directed or undirected.  `f` should be a
% function that takes an adjacency matrix (and optionally a parameter
% `a`) as its input(s) and returns a 1xN community affiliation matrix
% (i.e. a vector where each node is associated with an identifier, and
% nodes with the same identifier are in the same community).  The next
% optional argument is `it`, which is the number of iterations over
% which to calculate the versatility; usually the defualt of 1000 is
% more than enough, but this number can be reduced to around 100 for
% good approximations.  The last optional argument is `processors`,
% which should be greater than 1 if you would like to use
% parallelization, or 1 if you either a) use Octave or b) don't
% have the parallelization toolbox.
%
% This function will return the versatility of each node in the
% network as a 1xN vector.
%
% As you likely know, versatility describes how closely each node is
% affiliated with a module.  It does not depend on any particular
% modular structure; it only depends on the algorithm and its
% parametrization.  Averaging a node's versatility across many
% parameters creates a stable measure that is consistent across
% different algorithms in most networks, and also eliminates the
% dependence on the parameterization.
%
% This default argument for this function comes from the Brain
% Connectivity Toolbox for the "community_louvain" function, but it
% does not necessarily depend on this package.
%
% Example usage:
%
%     % Louvain versatility at gamma=2.0
%     find_nodal_versatility(G, @community_louvain, 2.0)
%     % Louvain versatility at gamma=1.0, run on 3 processors
%     find_nodal_versatility(G, [], [], [], 3)
%
% Copyright 2016 Max Shinn <maxwell.shinn@yale.edu>
% Available under the GPLv3

% Default arguments: `g` (N/A), `f` (@community_louvain), `a` (no
% extra argument to `f`, or 1.0 if `f` is set to the default), `it`
% (1000).
function V = find_nodal_versatility(g, f, a, it, processors)
    % Set default arguments
    if ~exist('g','var') || isempty(g) || all(size(g) ~= size(g'))
        error('Invalid matrix');
    end
    if ~exist('f', 'var') || isempty(f)
        f = @community_louvain;
        % Only assign a=1.0 if we're using community_louvain
        if ~exist('a', 'var') || isempty(ar)
            a = 1.0;
        end
    end
    if ~exist('a', 'var') || isempty(a) % If `f` was assigned but `a` wasn't
        a = 'invalid'; % A hack to make things easier
    end
    if ~exist('it','var') || isempty(it)
        it = 1000;
    end
    if ~exist('processors','var') || isempty(processors)
        processors = 1;
    end
    
    if processors == 1
        CM = consensus_matrix(g, f, a, it);
    else
        CM = consensus_matrix_par(g, f, a, it, processors);
    end
    
    Cs = sin(pi*CM);
    V = sum(Cs, 1)./sum(CM, 1); % CM/Cs are symmetric so axis doesn't matter
    V(V<1e-10) = 0; % Stupid floats
end

% consensus_matrix finds the consensus matrix (i.e. the association
% matrix) of the matrix `g`.  The matrix M returned by this function
% satisfies the property such that, for two nodes i and j, M_{i,j} is
% the probability that i and j will be in the same module on a random
% run of the function `f`.  The arguments should satisfy the same
% properties as in find_nodal_versatility with one exception: `a`
% should never be empty.  Instead, if no argument is to be passed to
% the function `f`, it should be set to the string 'invalid'.
function M = consensus_matrix(g, f, a, it)
    sg = size(g);
    M = zeros(sg);
    if a == 'invalid'
        for i=[1:it]
            p = f(g); % Get community affiliations
            p = p(:); % Assert column vector
            M = M + (repmat(p, 1, sg(1)) == repmat(p', sg(2), 1)); % Matrix where M_{i,j} == 1 iff p_i == p_j
        end
    else
        for i=[1:it]
            p = f(g, a); % Get community affiliations
            p = p(:); % Assert column vector
            M = M + (repmat(p, 1, sg(1)) == repmat(p', sg(2), 1)); % Matrix where M_{i,j} == 1 iff p_i == p_j
        end
    end
    M = M / it;
end

% Almost exactly the same as consensus_matrix.  However, it will
% run the process in parallel.  It uses more memory than
% consensus_matrix, but is usually faster.  The extra argument
% `processors` specifies the number of processors to use in the
% parallel pool.  It should always be an integer greater than 1.
function M = consensus_matrix_par(g, f, a, it, processors)
    if exist('parpool') == 0
        error('No parallelization support')
    end
    % Create parallel pool if it doesn't exist
    p = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(p)
        parpool(processors);
    elseif p.NumWorkers ~= processors % If existing pool is the wrong size, make a new one
        delete(p);
        parpool(processors);
    end
    
    sg = size(g);
    cons = zeros(sg(1), sg(2), processors);
    if a == 'invalid'
        parfor i=[1:processors]
            conssub = zeros(sg(1), sg(2))
            for j=[1:ceil(it/processors)]
                p = f(g); % Get community affiliations
                p = p(:); % Assert column vector
                conssub = conssub + (repmat(p, 1, sg(1)) == repmat(p', sg(2), 1)); % Matrix where M_{i,j} == 1 iff p_i == p_j
            end
            cons(:,:,i) = conssub/ceil(it/processors)
        end
    else
        parfor i=[1:processors]
            conssub = zeros(sg(1), sg(2))
            for j=[1:ceil(it/processors)]
                p = f(g, a); % Get community affiliations
                p = p(:); % Assert column vector
                conssub = conssub + (repmat(p, 1, sg(1)) == repmat(p', sg(2), 1)); % Matrix where M_{i,j} == 1 iff p_i == p_j
            end
            cons(:,:,i) = conssub/ceil(it/processors)
        end
    end
    M = mean(cons, 3);
end

