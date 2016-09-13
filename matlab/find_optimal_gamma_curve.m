% find_optimal_gamma_curve will plot versatility across different
% gamma values using a specified community detection algorithm.  `g`
% is an NxN adjacency matrix, which can either be weighted or
% unweighted, directed or undirected.  `f` should be a function that
% takes an adjacency matrix as its first argument and a resolution
% parameter as its input(s) and returns a 1xN community affiliation
% matrix (i.e. a vector where each node is associated with an
% identifier, and nodes with the same identifier are in the same
% community).  The argument is `ar`, which should be the values of the
% resolution parameter to plot.  The next optional argument is `it`,
% which is the number of iterations over which to calculate the
% versatility AT EACH GAMMA; the defualt of 100 gives about .01
% variance on average.  The last optional argument is `processors`,
% which should be 1 if you would like to use parallelization.
%
% This returns the computed versatility values, corresponding to
% each value of `ar`.
% 
% Example usage:
%
%     % Display the curve from .1 to 4.0
%     find_optimal_gamma_curve(G)
%     % Display the curve from .1 to 4.0, run on 3 processors
%     find_optimal_gamma_curve(G, [], [], [], 3)
%     % Display the curve with a different algorithm from .5 to 2.0, using 200 iterations per gamma
%     find_optimal_gamma_curve(G, @modularity_louvain_dir, .5:.05:2.0, 200)
%
%
% Copyright 2016 Max Shinn <mws41@cam.ac.uk>
% Available under the GPLv3

% Default arguments: `g` (N/A), `f` (@community_louvain), `a` (.1 to
% 4.0, spaced by .1 intervals), `it` (100), `processors` (1)
function vs = find_optimal_gamma_curve(g, f, ar, it, processors)
    % Set default arguments
    if ~exist('g','var') || isempty(g) || all(size(g) ~= size(g'))
        error('Invalid matrix');
    end
    if ~exist('f', 'var') || isempty(f)
        f = @community_louvain;
        % Only assign a=1.0 if we're using community_louvain
        if ~exist('ar', 'var') || isempty(ar)
            ar = .1:.1:4;
        end
    end
    if ~exist('ar', 'var') || isempty(ar) % If `f` was assigned but `a` wasn't
        error('The argument `ar` must be assigned if a function is given');
    end
    if ~exist('it','var') || isempty(it)
        it = 100;
    end
    if ~exist('processors','var') || isempty(processors)
        processors = 1;
    end

    % Do the grunt work
    vs = [];
    for a=ar
        v = mean(find_nodal_versatility(g, f, a, it, processors));
        vs = [vs v];
        disp(['Processed ', num2str(a)])
    end

    plot(ar, vs)
    xlabel('Resolution parameters')
    ylabel('Versatility')
    title('Parameters that minimize versatility make the best communities')
    
end