function out = getBSplinePP( order )
    % x in range [a, b] -- BIG ASSUMPTION -- will return non-zero values
    % outside of [a,b]
    % create bspline on "normal" interval
    switch order
        case (1)
            out = mkpp([-0.5, 0.5], [1]);
        case (2)
            out = mkpp([-1, 0, 1], [[1,0];[-1,1]]);
        case (3)
            brk = -1.5:1:1.5;
            out = mkpp(brk, [[0.5, 0, 0];[-1, 1, 0.5];[0.5, -1, 0.5]]);
        case (4)
            brk = -2:1:2;
            out = mkpp(brk, [[1/6.0, 0, 0, 0];[-0.5, 0.5, 0.5, 1.0/6.0];[0.5, -1, 0, 2/3];[-1/6.0, 0.5, -0.5, 1.0/6.0]] );
        case (5)
            brk = -2.5:2.5;
            out = mkpp(brk, [[1 / 24.0, 0, 0, 0, 0];
                    [-1/6, 1/6, 0.25, 1/6, 1/24];
                    [0.25, -0.5, -0.25, 0.5, 11 / 24];
                    [-1/6, 0.5, -0.25, -0.5, 11/ 24];
                    [1/ 24, -1/6.0, 0.25, -1/6, 1/24]] );
        case (6)
            brk = -3:3;
            out = mkpp(brk, [[1/120, 0, 0, 0, 0, 0];
                    [-1/24, 1/24, 1/12, 1/12, 1/24, 1/120];
                    [1/12, -1/6, -1/6, 1/6, 5/12, 13/60];
                    [-1/12, 0.25, 0, -1/2, 0, 11/20];
                    [1/24, -1/6, 1/6, 1/6, -5/12, 13/60];
                    [-1/120, 1/24, -1/12, 1/12, -1/24, 1/120]] );
        case (7)
            brk = -3.5:3.5;
            out = mkpp(brk, [[1/720, 0, 0, 0, 0, 0, 0];
                    [-1/120, 1/120, 1/48, 1/36, 1/48, 1/120, 1/720 ];
                    [1/48, -1/24, -1/16, 1/36, 3/16, 5/24, 19/240];
                    [-1/36, 1/12, 1/24, -2/9, -5/24, 1/3, 151/360];
                    [1/48, -1/12, 1/24, 2/9, -5/24, -1/3, 151/360];
                    [-1/120, 1/24, -1/16, -1/36, 3/16, -5/24, 19/240];
                    [1/720, -1/120, 1/48, -1/36, 1/48, -1/120, 1/720]] );
        case(8)
            brk = -4:4;
            out = mkpp( brk, [[1/5040, 0, 0, 0, 0, 0, 0, 0];
                              [-1/720, 1/720, 1/240, 1/144, 1/144, 1/240, 1/720, 1/5040 ];
                              [1/240, -1/120, -1/60, 0, 1/18, 1/10, 7/90, 1/42];
                              [-1/144, 1/48, 1/48, -1/16, -19/144, 1/16, 49/144, 397/1680];
                              [1/144, -1/36, 0, 1/9, 0, -1/3, 0, 151/315 ];
                              [-1/240, 1/48, -1/48, -1/16, 19/144, 1/16, -49/144, 397/1680];
                              [ 1/720, -1/120, 1/60, 0, -1/18, 1/10, -7/90, 1/42 ];
                              [ -1/5040, 1/720, -1/240, 1/144, -1/144, 1/240, -1/720, 1/5040]] );
         case(9)
            brk = -4.5:4.5;
            out = mkpp( brk, [[1/40320, 0, 0, 0, 0, 0, 0, 0, 0];
                              [-1/5040, -1/5040, 1/1440, 1/720, 1/576, 1/720, 1/1440, 1/5040, 1/40320];
                              [1/1440, -1/720, -1/288, -1/720, 7/576, 23/720, 11/288, 17/720, 247/40320];
                              [-1/720, 1/240, 1/160, -1/80, -3/64, -1/80, 21/160, 17/80, 477/4480];
                              [1/576, -1/144, -1/288, 5/144, 19/576, -19/144, -49/288, 35/144, 15619/40320];
                              [ -1/720, 1/144, -1/288, -5/144, 19/576, 19/144, -49/288, -35/144, 15619/40320];
                              [ 1/1440, -1/240, 1/160, 1/80, -3/64, 1/80, 21/160, -17/80 477/4480];
                              [ -1/5040, 1/720, -1/288, 1/720, 7/576, -23/720, 11/288, -17/720, 247/40320 ];
                              [ 1/40320, -1/5040, 1/1440, -1/720, 1/576, -1/720, 1/1440, -1/5040, 1/40320]] );
         case(10)
            brk = -5:5;
            out = mkpp( brk, [[1/362880, 0, 0, 0, 0, 0, 0, 0, 0,0];
                              [-1/40320, 1/40320, 1/10080, 1/4320, 1/2880, 1/2880, 1/4320, 1/10080, 1/40320, 1/362880];
                              [1/10080, -1/5040, -1/1680, -1/2160, 1/480, 11/1440, 1/80, 59/5040, 41/6720, 251/181440];
                              [-1/4320, 1/1440, 1/720, -1/540, -17/1440, -1/90, 67/2160, 17/180, 289/2880, 913/22680];
                              [1/2880, -1/720, -1/720, 17/2160, 23/1440, -43/1440, -217/2160, 11/720, 809/2880, 44117/181440];
                              [-1/2880, 1/576, 0, -5/432, 0, 19/288, 0, -35/144, 0, 15619/36288];
                              [1/4320, -1/720, 1/720, 17/2160, -23/1440, -43/1440, 217/2160, 11/720, -809/2880, 44117/181440];
                              [-1/10080, 1/1440, -1/720, -1/540, 17/1440, -1/90, -67/2160, 17/180, -289/2880, 913/22680];
                              [1/40320, -1/5040, 1/1680, -1/2160, -1/480, 11/1440, -1/80, 59/5040, -41/6720, 251/181440];
                              [-1/362880, 1/40320, -1/10080, 1/4320, -1/2880, 1/2880, -1/4320, 1/10080, -1/40320, 1/362880]] );
                          
        otherwise
            fprintf('Order %d is not supported\n', order);
    end 