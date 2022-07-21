ifort -shared -fPIC -fpp -extend-source 132 -threads -static-intel -O3 -o GeoMaterials64.so \
../code/String.f90            \
../code/GlobalConstants.f90   \
../code/Debug.f90             \
../code/MathLib.f90           \
../code/StressStrainLib.f90   \
../code/YieldSurfaceLib.f90   \
../code/SolversLib.f90        \
../code/ElastoPlastic.f90     \
../code/ViscoPlastic.f90      \
../code/userAdd.f90           \
../code/userMod.f90

