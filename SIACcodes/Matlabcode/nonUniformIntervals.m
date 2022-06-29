function bpoints = nonUniformIntervals( numIntervals, domain, randVal )
    %NONUNIFORMINTERVALS    Creates a vector of breakpoings for specified number of intervals over specified domain
    %   NONUNIFORMINTERVALS( N, D, R ) is a vector with N + 1 breakpoints.
    %   They occur in the domain [0, D].  If R is true, then the values are
    %   randomly distributed in the domain.  If false, it always returns
    %   the same vector
    
    if ( ~randVal )
        fprintf('Resetting random\n');
        rand('state', 0 );
    end
    
    last = numIntervals + 1;
    bpoints( 1, 2:last ) = rand(1, numIntervals );
    for i = 2:last
        bpoints(i) = bpoints(i-1) + bpoints(i);
    end
    scale = domain / bpoints( last );
    bpoints = bpoints * scale;
    