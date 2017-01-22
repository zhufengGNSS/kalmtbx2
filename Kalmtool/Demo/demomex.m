%---------- PC/Linux ----------
if strcmpi(computer,'LNX86'),
   mex dd1fall.c ../kalmlblx.o
   mex dd2fall.c ../kalmlblx.o
   mex dd1agv.c ../kalmlblx.o
   mex dd1magv.c ../kalmlblx.o
   mex dd2agv.c ../kalmlblx.o
   mex dd2magv.c ../kalmlblx.o

%---------- MS-Windows ----------
elseif strcmpi(computer,'PCWIN'),
        mex dd1fall.c ../kalmlblcc.obj
        mex dd2fall.c ../kalmlblcc.obj
        mex dd1agv.c ../kalmlblcc.obj
        mex dd1magv.c ../kalmlblcc.obj
        mex dd2agv.c ../kalmlblcc.obj
        mex dd2magv.c ../kalmlblcc.obj
end
