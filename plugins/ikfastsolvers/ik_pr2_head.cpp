#define IKFAST_NAMESPACE ik_pr2_head
#include "plugindefs.h"

/// autogenerated analytical inverse kinematics code from ikfast program part of OpenRAVE
/// \author Rosen Diankov
///
/// Licensed under the Apache License, Version 2.0 (the "License");
/// you may not use this file except in compliance with the License.
/// You may obtain a copy of the License at
///     http://www.apache.org/licenses/LICENSE-2.0
/// 
/// Unless required by applicable law or agreed to in writing, software
/// distributed under the License is distributed on an "AS IS" BASIS,
/// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
/// See the License for the specific language governing permissions and
/// limitations under the License.
///
/// ikfast version 48 generated on 2011-10-14 18:28:26.655655
/// To compile with gcc:
///     gcc -lstdc++ ik.cpp
/// To compile without any main function as a shared object:
///     gcc -fPIC -lstdc++ -DIKFAST_NO_MAIN -shared -Wl,-soname,ik.so -o ik.so ik.cpp
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>
#include <complex>

#ifdef IKFAST_HEADER
#include IKFAST_HEADER
#endif

#ifndef IKFAST_ASSERT
#include <stdexcept>
#include <sstream>

#ifdef _MSC_VER
#ifndef __PRETTY_FUNCTION__
#define __PRETTY_FUNCTION__ __FUNCDNAME__
#endif
#endif

#ifndef __PRETTY_FUNCTION__
#define __PRETTY_FUNCTION__ __func__
#endif

#define IKFAST_ASSERT(b) { if( !(b) ) { std::stringstream ss; ss << "ikfast exception: " << __FILE__ << ":" << __LINE__ << ": " <<__PRETTY_FUNCTION__ << ": Assertion '" << #b << "' failed"; throw std::runtime_error(ss.str()); } }

#endif

#if defined(_MSC_VER)
#define IKFAST_ALIGNED16(x) __declspec(align(16)) x
#else
#define IKFAST_ALIGNED16(x) x __attribute((aligned(16)))
#endif

#define IK2PI  ((IKReal)6.28318530717959)
#define IKPI  ((IKReal)3.14159265358979)
#define IKPI_2  ((IKReal)1.57079632679490)

#ifdef _MSC_VER
#ifndef isnan
#define isnan _isnan
#endif
#endif // _MSC_VER

// defined when creating a shared object/dll
#ifdef IKFAST_CLIBRARY
#ifdef _MSC_VER
#define IKFAST_API extern "C" __declspec(dllexport)
#else
#define IKFAST_API extern "C"
#endif
#else
#define IKFAST_API
#endif

// lapack routines
extern "C" {
  void dgetrf_ (const int* m, const int* n, double* a, const int* lda, int* ipiv, int* info);
  void zgetrf_ (const int* m, const int* n, std::complex<double>* a, const int* lda, int* ipiv, int* info);
  void dgetri_(const int* n, const double* a, const int* lda, int* ipiv, double* work, const int* lwork, int* info);
  void dgesv_ (const int* n, const int* nrhs, double* a, const int* lda, int* ipiv, double* b, const int* ldb, int* info);
  void dgetrs_(const char *trans, const int *n, const int *nrhs, double *a, const int *lda, int *ipiv, double *b, const int *ldb, int *info);
  void dgeev_(const char *jobvl, const char *jobvr, const int *n, double *a, const int *lda, double *wr, double *wi,double *vl, const int *ldvl, double *vr, const int *ldvr, double *work, const int *lwork, int *info);
}

using namespace std; // necessary to get std math routines

#ifdef IKFAST_NAMESPACE
namespace IKFAST_NAMESPACE {
#endif

#ifdef IKFAST_REAL
typedef IKFAST_REAL IKReal;
#else
typedef double IKReal;
#endif

class IKSolution
{
public:
    /// Gets a solution given its free parameters
    /// \param pfree The free parameters required, range is in [-pi,pi]
    void GetSolution(IKReal* psolution, const IKReal* pfree) const {
        for(std::size_t i = 0; i < basesol.size(); ++i) {
            if( basesol[i].freeind < 0 )
                psolution[i] = basesol[i].foffset;
            else {
                IKFAST_ASSERT(pfree != NULL);
                psolution[i] = pfree[basesol[i].freeind]*basesol[i].fmul + basesol[i].foffset;
                if( psolution[i] > IKPI ) {
                    psolution[i] -= IK2PI;
                }
                else if( psolution[i] < -IKPI ) {
                    psolution[i] += IK2PI;
                }
            }
        }
    }

    /// Gets the free parameters the solution requires to be set before a full solution can be returned
    /// \return vector of indices indicating the free parameters
    const std::vector<int>& GetFree() const { return vfree; }

    struct VARIABLE
    {
        VARIABLE() : freeind(-1), fmul(0), foffset(0) {}
        VARIABLE(int freeind, IKReal fmul, IKReal foffset) : freeind(freeind), fmul(fmul), foffset(foffset) {}
        int freeind;
        IKReal fmul, foffset; ///< joint value is fmul*sol[freeind]+foffset
    };

    std::vector<VARIABLE> basesol;       ///< solution and their offsets if joints are mimiced
    std::vector<int> vfree;
};

inline float IKabs(float f) { return fabsf(f); }
inline double IKabs(double f) { return fabs(f); }

inline float IKlog(float f) { return logf(f); }
inline double IKlog(double f) { return log(f); }

#ifndef IKFAST_SINCOS_THRESH
#define IKFAST_SINCOS_THRESH ((IKReal)0.000001)
#endif

#ifndef IKFAST_ATAN2_MAGTHRESH
#define IKFAST_ATAN2_MAGTHRESH ((IKReal)2e-6)
#endif

inline float IKasin(float f)
{
IKFAST_ASSERT( f > -1-IKFAST_SINCOS_THRESH && f < 1+IKFAST_SINCOS_THRESH ); // any more error implies something is wrong with the solver
if( f <= -1 ) return -IKPI_2;
else if( f >= 1 ) return IKPI_2;
return asinf(f);
}
inline double IKasin(double f)
{
IKFAST_ASSERT( f > -1-IKFAST_SINCOS_THRESH && f < 1+IKFAST_SINCOS_THRESH ); // any more error implies something is wrong with the solver
if( f <= -1 ) return -IKPI_2;
else if( f >= 1 ) return IKPI_2;
return asin(f);
}

// return positive value in [0,y)
inline float IKfmod(float x, float y)
{
    while(x < 0) {
        x += y;
    }
    return fmodf(x,y);
}

// return positive value in [0,y)
inline float IKfmod(double x, double y)
{
    while(x < 0) {
        x += y;
    }
    return fmod(x,y);
}

inline float IKacos(float f)
{
IKFAST_ASSERT( f > -1-IKFAST_SINCOS_THRESH && f < 1+IKFAST_SINCOS_THRESH ); // any more error implies something is wrong with the solver
if( f <= -1 ) return IKPI;
else if( f >= 1 ) return 0;
return acosf(f);
}
inline double IKacos(double f)
{
IKFAST_ASSERT( f > -1-IKFAST_SINCOS_THRESH && f < 1+IKFAST_SINCOS_THRESH ); // any more error implies something is wrong with the solver
if( f <= -1 ) return IKPI;
else if( f >= 1 ) return 0;
return acos(f);
}
inline float IKsin(float f) { return sinf(f); }
inline double IKsin(double f) { return sin(f); }
inline float IKcos(float f) { return cosf(f); }
inline double IKcos(double f) { return cos(f); }
inline float IKtan(float f) { return tanf(f); }
inline double IKtan(double f) { return tan(f); }
inline float IKsqrt(float f) { if( f <= 0.0f ) return 0.0f; return sqrtf(f); }
inline double IKsqrt(double f) { if( f <= 0.0 ) return 0.0; return sqrt(f); }
inline float IKatan2(float fy, float fx) {
    if( isnan(fy) ) {
        IKFAST_ASSERT(!isnan(fx)); // if both are nan, probably wrong value will be returned
        return IKPI_2;
    }
    else if( isnan(fx) ) {
        return 0;
    }
    return atan2f(fy,fx);
}
inline double IKatan2(double fy, double fx) {
    if( isnan(fy) ) {
        IKFAST_ASSERT(!isnan(fx)); // if both are nan, probably wrong value will be returned
        return IKPI_2;
    }
    else if( isnan(fx) ) {
        return 0;
    }
    return atan2(fy,fx);
}

inline float IKsign(float f) {
    if( f > 0 ) {
        return 1.0f;
    }
    else if( f < 0 ) {
        return -1.0f;
    }
    return 0;
}

inline double IKsign(double f) {
    if( f > 0 ) {
        return 1.0;
    }
    else if( f < 0 ) {
        return -1.0;
    }
    return 0;
}

/// solves the inverse kinematics equations.
/// \param pfree is an array specifying the free joints of the chain.
IKFAST_API void fk(const IKReal* j, IKReal* eetrans, IKReal* eerot) {
IKReal x0,x1,x2,x3;
x0=IKcos(j[0]);
x1=IKcos(j[1]);
x2=IKsin(j[0]);
x3=IKsin(j[1]);
eetrans[0]=((-0.0170700000000000)+(((0.0232000000000000)*(x0)*(x1)))+(((-0.0300000000000000)*(x2)))+(((0.0680000000000000)*(x0)))+(((0.0980000000000000)*(x0)*(x3))));
eetrans[1]=((((0.0980000000000000)*(x2)*(x3)))+(((0.0300000000000000)*(x0)))+(((0.0680000000000000)*(x2)))+(((0.0232000000000000)*(x1)*(x2))));
eetrans[2]=((0.381450000000000)+(((0.0980000000000000)*(x1)))+(((-0.0232000000000000)*(x3))));
eerot[0]=((x0)*(x1));
eerot[1]=((x1)*(x2));
eerot[2]=((-1.00000000000000)*(x3));
}

IKFAST_API int getNumFreeParameters() { return 0; }
IKFAST_API int* getFreeParameters() { return NULL; }
IKFAST_API int getNumJoints() { return 2; }

IKFAST_API int getIKRealSize() { return sizeof(IKReal); }

IKFAST_API int getIKType() { return 0x23000006; }

class IKSolver {
public:
IKReal j13,cj13,sj13,htj13,j14,cj14,sj14,htj14,new_px,px,npx,new_py,py,npy,new_pz,pz,npz,pp;

bool ik(const IKReal* eetrans, const IKReal* eerot, const IKReal* pfree, std::vector<IKSolution>& vsolutions) {
for(int dummyiter = 0; dummyiter < 1; ++dummyiter) {
vsolutions.resize(0); vsolutions.reserve(8);
px = eetrans[0]; py = eetrans[1]; pz = eetrans[2];

new_px=((0.0170700000000000)+(px));
new_py=py;
new_pz=((-0.381450000000000)+(pz));
px = new_px; py = new_py; pz = new_pz;
pp=(((px)*(px))+((py)*(py))+((pz)*(pz)));
{
IKReal dummyeval[1];
dummyeval[0]=((1.00000000000000)+(((33.3333333333333)*(py))));
if( IKabs(dummyeval[0]) < 0.0000010000000000  )
{
continue;

} else
{
{
IKReal j13array[2], cj13array[2], sj13array[2];
bool j13valid[2]={false};
IKReal x4=((2.00000000000000)*(py));
IKReal x5=((0.0600000000000000)+(x4));
IKReal x6=((IKabs(x5) != 0)?((IKReal)1/(x5)):(IKReal)1.0e30);
IKReal x7=(px)*(px);
IKReal x8=((16.0000000000000)*(x7));
IKReal x9=(py)*(py);
IKReal x10=((16.0000000000000)*(x9));
IKReal x11=((-0.0144000000000000)+(x8)+(x10));
if( (x11) < (IKReal)-0.00001 )
    continue;
IKReal x12=IKsqrt(x11);
IKReal x13=((0.500000000000000)*(x12)*(x6));
IKReal x14=((2.00000000000000)*(px)*(x6));
j13array[0]=((2.00000000000000)*(atan(((((-1.00000000000000)*(x14)))+(x13)))));
sj13array[0]=IKsin(j13array[0]);
cj13array[0]=IKcos(j13array[0]);
j13array[1]=((-2.00000000000000)*(atan(((x13)+(x14)))));
sj13array[1]=IKsin(j13array[1]);
cj13array[1]=IKcos(j13array[1]);
if( j13array[0] > IKPI )
{
    j13array[0]-=IK2PI;
}
else if( j13array[0] < -IKPI )
{    j13array[0]+=IK2PI;
}
j13valid[0] = true;
if( j13array[1] > IKPI )
{
    j13array[1]-=IK2PI;
}
else if( j13array[1] < -IKPI )
{    j13array[1]+=IK2PI;
}
j13valid[1] = true;
if( j13valid[0] && j13valid[1] && IKabs(cj13array[0]-cj13array[1]) < 0.0001 && IKabs(sj13array[0]-sj13array[1]) < 0.0001 )
{
    j13valid[1]=false;
}
for(int ij13 = 0; ij13 < 2; ++ij13)
{
if( !j13valid[ij13] )
{
    continue;
}
j13 = j13array[ij13]; cj13 = cj13array[ij13]; sj13 = sj13array[ij13];

{
IKReal dummyeval[1];
IKReal x15=(sj13)*(sj13);
dummyeval[0]=((((5.13777777777778)*(x15)))+((cj13)*(cj13))+(((1111.11111111111)*((py)*(py))))+(((1111.11111111111)*(x15)*((pz)*(pz))))+(((4.53333333333333)*(cj13)*(sj13)))+(((-151.111111111111)*(py)*(sj13)))+(((-66.6666666666667)*(cj13)*(py))));
if( IKabs(dummyeval[0]) < 0.0000010000000000  )
{
{
IKReal dummyeval[1];
IKReal x16=(cj13)*(cj13);
IKReal x17=(pz)*(pz);
dummyeval[0]=((1.00000000000000)+(((-1.00000000000000)*(x16)))+(((104.123281965848)*(x16)*(x17)))+(((-104.123281965848)*(x17))));
if( IKabs(dummyeval[0]) < 0.0000001000000000  )
{
{
IKReal dummyeval[1];
IKReal x18=(cj13)*(cj13);
IKReal x19=(px)*(px);
dummyeval[0]=((1.00000000000000)+(((1111.11111111111)*(x18)*(x19)))+(((-1111.11111111111)*(x19)))+(((2222.22222222222)*(cj13)*(px)*(py)*(sj13)))+(((-1111.11111111111)*(x18)*((py)*(py)))));
if( IKabs(dummyeval[0]) < 0.0000010000000000  )
{
{
IKReal evalcond[1];
evalcond[0]=((-3.14159265358979)+(IKfmod(((3.14159265358979)+(j13)), 6.28318530717959)));
if( IKabs(evalcond[0]) < 0.0000010000000000  )
{
{
IKReal dummyeval[1];
dummyeval[0]=((1.00000000000000)+(((-1111.11111111111)*((py)*(py)))));
if( IKabs(dummyeval[0]) < 0.0000010000000000  )
{
continue;

} else
{
{
IKReal j14array[4], cj14array[4], sj14array[4];
bool j14valid[4]={false};
IKReal x20=(py)*(py);
IKReal x21=((0.000900000000000000)+(((-1.00000000000000)*(x20))));
IKReal x22=((IKabs(x21) != 0)?((IKReal)1/(x21)):(IKReal)1.0e30);
IKReal x23=((IKabs(x21) != 0)?(pow(x21,-0.500000000000000)):(IKReal)1.0e30);
if( (x21) < (IKReal)-0.00001 )
    continue;
IKReal x24=IKsqrt(x21);
cj14array[0]=((-1.00000000000000)*(x23)*(x24));
cj14array[2]=((x23)*(x24));
if( cj14array[0] >= -1-IKFAST_SINCOS_THRESH && cj14array[0] <= 1+IKFAST_SINCOS_THRESH )
{
    j14valid[0] = j14valid[1] = true;
    j14array[0] = IKacos(cj14array[0]);
    sj14array[0] = IKsin(j14array[0]);
    cj14array[1] = cj14array[0];
    j14array[1] = -j14array[0];
    sj14array[1] = -sj14array[0];
}
else if( isnan(cj14array[0]) )
{
    // probably any value will work
    j14valid[0] = true;
    cj14array[0] = 1; sj14array[0] = 0; j14array[0] = 0;
}
if( cj14array[2] >= -1-IKFAST_SINCOS_THRESH && cj14array[2] <= 1+IKFAST_SINCOS_THRESH )
{
    j14valid[2] = j14valid[3] = true;
    j14array[2] = IKacos(cj14array[2]);
    sj14array[2] = IKsin(j14array[2]);
    cj14array[3] = cj14array[2];
    j14array[3] = -j14array[2];
    sj14array[3] = -sj14array[2];
}
else if( isnan(cj14array[2]) )
{
    // probably any value will work
    j14valid[2] = true;
    cj14array[2] = 1; sj14array[2] = 0; j14array[2] = 0;
}
if( j14valid[0] && j14valid[1] && IKabs(cj14array[0]-cj14array[1]) < 0.0001 && IKabs(sj14array[0]-sj14array[1]) < 0.0001 )
{
    j14valid[1]=false;
}
if( j14valid[0] && j14valid[2] && IKabs(cj14array[0]-cj14array[2]) < 0.0001 && IKabs(sj14array[0]-sj14array[2]) < 0.0001 )
{
    j14valid[2]=false;
}
if( j14valid[0] && j14valid[3] && IKabs(cj14array[0]-cj14array[3]) < 0.0001 && IKabs(sj14array[0]-sj14array[3]) < 0.0001 )
{
    j14valid[3]=false;
}
if( j14valid[1] && j14valid[2] && IKabs(cj14array[1]-cj14array[2]) < 0.0001 && IKabs(sj14array[1]-sj14array[2]) < 0.0001 )
{
    j14valid[2]=false;
}
if( j14valid[1] && j14valid[3] && IKabs(cj14array[1]-cj14array[3]) < 0.0001 && IKabs(sj14array[1]-sj14array[3]) < 0.0001 )
{
    j14valid[3]=false;
}
if( j14valid[2] && j14valid[3] && IKabs(cj14array[2]-cj14array[3]) < 0.0001 && IKabs(sj14array[2]-sj14array[3]) < 0.0001 )
{
    j14valid[3]=false;
}
for(int ij14 = 0; ij14 < 4; ++ij14)
{
if( !j14valid[ij14] )
{
    continue;
}
j14 = j14array[ij14]; cj14 = cj14array[ij14]; sj14 = sj14array[ij14];

IKReal soleval[1];
soleval[0]=((-0.0232000000000000)+(((-1.00000000000000)*(cj14)*(((0.0680000000000000)+(((-1.00000000000000)*(cj13)*(px)))+(((-1.00000000000000)*(py)*(sj13)))))))+(((-1.00000000000000)*(pz)*(sj14))));
if( soleval[0] > 0.0000000000000000  )
{
vsolutions.push_back(IKSolution()); IKSolution& solution = vsolutions.back();
solution.basesol.resize(2);
solution.basesol[0].foffset = j13;
solution.basesol[1].foffset = j14;
solution.vfree.resize(0);
}
}
}

}

}

} else
{
evalcond[0]=((-3.14159265358979)+(IKfmod(j13, 6.28318530717959)));
if( IKabs(evalcond[0]) < 0.0000010000000000  )
{
{
IKReal dummyeval[1];
dummyeval[0]=((1.00000000000000)+(((-1111.11111111111)*((py)*(py)))));
if( IKabs(dummyeval[0]) < 0.0000010000000000  )
{
continue;

} else
{
{
IKReal j14array[4], cj14array[4], sj14array[4];
bool j14valid[4]={false};
IKReal x25=(py)*(py);
IKReal x26=((0.000900000000000000)+(((-1.00000000000000)*(x25))));
IKReal x27=((IKabs(x26) != 0)?((IKReal)1/(x26)):(IKReal)1.0e30);
IKReal x28=((IKabs(x26) != 0)?(pow(x26,-0.500000000000000)):(IKReal)1.0e30);
if( (x26) < (IKReal)-0.00001 )
    continue;
IKReal x29=IKsqrt(x26);
cj14array[0]=((-1.00000000000000)*(x28)*(x29));
cj14array[2]=((x28)*(x29));
if( cj14array[0] >= -1-IKFAST_SINCOS_THRESH && cj14array[0] <= 1+IKFAST_SINCOS_THRESH )
{
    j14valid[0] = j14valid[1] = true;
    j14array[0] = IKacos(cj14array[0]);
    sj14array[0] = IKsin(j14array[0]);
    cj14array[1] = cj14array[0];
    j14array[1] = -j14array[0];
    sj14array[1] = -sj14array[0];
}
else if( isnan(cj14array[0]) )
{
    // probably any value will work
    j14valid[0] = true;
    cj14array[0] = 1; sj14array[0] = 0; j14array[0] = 0;
}
if( cj14array[2] >= -1-IKFAST_SINCOS_THRESH && cj14array[2] <= 1+IKFAST_SINCOS_THRESH )
{
    j14valid[2] = j14valid[3] = true;
    j14array[2] = IKacos(cj14array[2]);
    sj14array[2] = IKsin(j14array[2]);
    cj14array[3] = cj14array[2];
    j14array[3] = -j14array[2];
    sj14array[3] = -sj14array[2];
}
else if( isnan(cj14array[2]) )
{
    // probably any value will work
    j14valid[2] = true;
    cj14array[2] = 1; sj14array[2] = 0; j14array[2] = 0;
}
if( j14valid[0] && j14valid[1] && IKabs(cj14array[0]-cj14array[1]) < 0.0001 && IKabs(sj14array[0]-sj14array[1]) < 0.0001 )
{
    j14valid[1]=false;
}
if( j14valid[0] && j14valid[2] && IKabs(cj14array[0]-cj14array[2]) < 0.0001 && IKabs(sj14array[0]-sj14array[2]) < 0.0001 )
{
    j14valid[2]=false;
}
if( j14valid[0] && j14valid[3] && IKabs(cj14array[0]-cj14array[3]) < 0.0001 && IKabs(sj14array[0]-sj14array[3]) < 0.0001 )
{
    j14valid[3]=false;
}
if( j14valid[1] && j14valid[2] && IKabs(cj14array[1]-cj14array[2]) < 0.0001 && IKabs(sj14array[1]-sj14array[2]) < 0.0001 )
{
    j14valid[2]=false;
}
if( j14valid[1] && j14valid[3] && IKabs(cj14array[1]-cj14array[3]) < 0.0001 && IKabs(sj14array[1]-sj14array[3]) < 0.0001 )
{
    j14valid[3]=false;
}
if( j14valid[2] && j14valid[3] && IKabs(cj14array[2]-cj14array[3]) < 0.0001 && IKabs(sj14array[2]-sj14array[3]) < 0.0001 )
{
    j14valid[3]=false;
}
for(int ij14 = 0; ij14 < 4; ++ij14)
{
if( !j14valid[ij14] )
{
    continue;
}
j14 = j14array[ij14]; cj14 = cj14array[ij14]; sj14 = sj14array[ij14];

IKReal soleval[1];
soleval[0]=((-0.0232000000000000)+(((-1.00000000000000)*(cj14)*(((0.0680000000000000)+(((-1.00000000000000)*(cj13)*(px)))+(((-1.00000000000000)*(py)*(sj13)))))))+(((-1.00000000000000)*(pz)*(sj14))));
if( soleval[0] > 0.0000000000000000  )
{
vsolutions.push_back(IKSolution()); IKSolution& solution = vsolutions.back();
solution.basesol.resize(2);
solution.basesol[0].foffset = j13;
solution.basesol[1].foffset = j14;
solution.vfree.resize(0);
}
}
}

}

}

} else
{
if( 1 )
{
continue;

} else
{
}
}
}
}

} else
{
{
IKReal j14array[4], cj14array[4], sj14array[4];
bool j14valid[4]={false};
IKReal x30=(px)*(px);
IKReal x31=(cj13)*(cj13);
IKReal x32=((2.00000000000000)*(cj13)*(px)*(py)*(sj13));
IKReal x33=((x30)*(x31));
IKReal x34=((0.000900000000000000)+(x33)+(x32));
IKReal x35=(py)*(py);
IKReal x36=((x31)*(x35));
IKReal x37=((x30)+(x36));
IKReal x38=((((-1.00000000000000)*(x37)))+(x34));
IKReal x39=((IKabs(x38) != 0)?((IKReal)1/(x38)):(IKReal)1.0e30);
IKReal x40=((IKabs(x38) != 0)?(pow(x38,-0.500000000000000)):(IKReal)1.0e30);
if( (x38) < (IKReal)-0.00001 )
    continue;
IKReal x41=IKsqrt(x38);
cj14array[0]=((-1.00000000000000)*(x40)*(x41));
cj14array[2]=((x40)*(x41));
if( cj14array[0] >= -1-IKFAST_SINCOS_THRESH && cj14array[0] <= 1+IKFAST_SINCOS_THRESH )
{
    j14valid[0] = j14valid[1] = true;
    j14array[0] = IKacos(cj14array[0]);
    sj14array[0] = IKsin(j14array[0]);
    cj14array[1] = cj14array[0];
    j14array[1] = -j14array[0];
    sj14array[1] = -sj14array[0];
}
else if( isnan(cj14array[0]) )
{
    // probably any value will work
    j14valid[0] = true;
    cj14array[0] = 1; sj14array[0] = 0; j14array[0] = 0;
}
if( cj14array[2] >= -1-IKFAST_SINCOS_THRESH && cj14array[2] <= 1+IKFAST_SINCOS_THRESH )
{
    j14valid[2] = j14valid[3] = true;
    j14array[2] = IKacos(cj14array[2]);
    sj14array[2] = IKsin(j14array[2]);
    cj14array[3] = cj14array[2];
    j14array[3] = -j14array[2];
    sj14array[3] = -sj14array[2];
}
else if( isnan(cj14array[2]) )
{
    // probably any value will work
    j14valid[2] = true;
    cj14array[2] = 1; sj14array[2] = 0; j14array[2] = 0;
}
if( j14valid[0] && j14valid[1] && IKabs(cj14array[0]-cj14array[1]) < 0.0001 && IKabs(sj14array[0]-sj14array[1]) < 0.0001 )
{
    j14valid[1]=false;
}
if( j14valid[0] && j14valid[2] && IKabs(cj14array[0]-cj14array[2]) < 0.0001 && IKabs(sj14array[0]-sj14array[2]) < 0.0001 )
{
    j14valid[2]=false;
}
if( j14valid[0] && j14valid[3] && IKabs(cj14array[0]-cj14array[3]) < 0.0001 && IKabs(sj14array[0]-sj14array[3]) < 0.0001 )
{
    j14valid[3]=false;
}
if( j14valid[1] && j14valid[2] && IKabs(cj14array[1]-cj14array[2]) < 0.0001 && IKabs(sj14array[1]-sj14array[2]) < 0.0001 )
{
    j14valid[2]=false;
}
if( j14valid[1] && j14valid[3] && IKabs(cj14array[1]-cj14array[3]) < 0.0001 && IKabs(sj14array[1]-sj14array[3]) < 0.0001 )
{
    j14valid[3]=false;
}
if( j14valid[2] && j14valid[3] && IKabs(cj14array[2]-cj14array[3]) < 0.0001 && IKabs(sj14array[2]-sj14array[3]) < 0.0001 )
{
    j14valid[3]=false;
}
for(int ij14 = 0; ij14 < 4; ++ij14)
{
if( !j14valid[ij14] )
{
    continue;
}
j14 = j14array[ij14]; cj14 = cj14array[ij14]; sj14 = sj14array[ij14];

IKReal soleval[1];
soleval[0]=((-0.0232000000000000)+(((-1.00000000000000)*(cj14)*(((0.0680000000000000)+(((-1.00000000000000)*(cj13)*(px)))+(((-1.00000000000000)*(py)*(sj13)))))))+(((-1.00000000000000)*(pz)*(sj14))));
if( soleval[0] > 0.0000000000000000  )
{
vsolutions.push_back(IKSolution()); IKSolution& solution = vsolutions.back();
solution.basesol.resize(2);
solution.basesol[0].foffset = j13;
solution.basesol[1].foffset = j14;
solution.vfree.resize(0);
}
}
}

}

}

} else
{
IKReal op[4+1], zeror[4];
int numroots;
IKReal x42=(cj13)*(cj13);
IKReal x43=(pz)*(pz);
IKReal x44=((0.0117600000000000)*(cj13)*(sj13));
IKReal x45=((0.0266560000000000)*(x42));
IKReal x46=((4.00000000000000)*(py)*(pz)*(sj13));
IKReal x47=((x42)*(x43));
IKReal x48=((0.00960400000000000)+(x47));
IKReal x49=((0.00960400000000000)*(x42));
IKReal x50=((x49)+(x43));
IKReal x51=((x48)+(((-1.00000000000000)*(x50))));
op[0]=x51;
op[1]=((0.0266560000000000)+(((-1.00000000000000)*(x45)))+(x46)+(x44));
op[2]=((0.0377040000000000)+(((-0.0341040000000000)*(x42)))+(((2.00000000000000)*(x43)))+(((-4.00000000000000)*((py)*(py))))+(((-2.00000000000000)*(x47)))+(((0.0163200000000000)*(cj13)*(sj13))));
op[3]=((0.0266560000000000)+(((-1.00000000000000)*(x45)))+(((-1.00000000000000)*(x46)))+(x44));
op[4]=x51;
polyroots4(op,zeror,numroots);
IKReal j14array[4], cj14array[4], sj14array[4], tempj14array[1];
int numsolutions = 0;
for(int ij14 = 0; ij14 < numroots; ++ij14)
{
IKReal htj14 = zeror[ij14];
tempj14array[0]=((2.00000000000000)*(atan(htj14)));
for(int kj14 = 0; kj14 < 1; ++kj14)
{
j14array[numsolutions] = tempj14array[kj14];
if( j14array[numsolutions] > IKPI )
{
    j14array[numsolutions]-=IK2PI;
}
else if( j14array[numsolutions] < -IKPI )
{
    j14array[numsolutions]+=IK2PI;
}
sj14array[numsolutions] = IKsin(j14array[numsolutions]);
cj14array[numsolutions] = IKcos(j14array[numsolutions]);
bool valid = true;
for( int kj14 = 0; kj14 < numsolutions; ++kj14)
{
    if( IKabs(cj14array[kj14]-cj14array[numsolutions]) < 0.0001 && IKabs(sj14array[kj14]-sj14array[numsolutions]) < 0.0001 )
    {
        valid=false; break;
    }
}
if( valid ) { numsolutions++; }
}
}
for(int ij14 = 0; ij14 < numsolutions; ++ij14)
    {
    j14 = j14array[ij14]; cj14 = cj14array[ij14]; sj14 = sj14array[ij14];

IKReal soleval[1];
soleval[0]=((-0.0232000000000000)+(((-1.00000000000000)*(cj14)*(((0.0680000000000000)+(((-1.00000000000000)*(cj13)*(px)))+(((-1.00000000000000)*(py)*(sj13)))))))+(((-1.00000000000000)*(pz)*(sj14))));
if( soleval[0] > 0.0000000000000000  )
{
vsolutions.push_back(IKSolution()); IKSolution& solution = vsolutions.back();
solution.basesol.resize(2);
solution.basesol[0].foffset = j13;
solution.basesol[1].foffset = j14;
solution.vfree.resize(0);
}
    }

}

}

} else
{
{
IKReal j14array[2], cj14array[2], sj14array[2];
bool j14valid[2]={false};
IKReal x52=((0.0680000000000000)*(sj13));
IKReal x53=((0.0300000000000000)*(cj13));
IKReal x54=((x53)+(x52));
IKReal x55=((py)+(((-1.00000000000000)*(x54))));
IKReal x56=(pz)*(pz);
IKReal x57=(sj13)*(sj13);
IKReal x58=((x56)*(x57));
IKReal x59=(x55)*(x55);
IKReal x60=((x59)+(x58));
if( (x60) < (IKReal)-0.00001 )
    continue;
IKReal x61=IKsqrt(x60);
IKReal x62=IKabs(x61);
IKReal x63=((IKabs(x62) != 0)?((IKReal)1/(x62)):(IKReal)1.0e30);
IKReal x64=((0.0980000000000000)*(sj13)*(x63));
if( (x64) < -1-IKFAST_SINCOS_THRESH || (x64) > 1+IKFAST_SINCOS_THRESH )
    continue;
IKReal x65=IKasin(x64);
IKReal x66=((pz)*(sj13));
if( IKabs(x66) < IKFAST_ATAN2_MAGTHRESH && IKabs(x55) < IKFAST_ATAN2_MAGTHRESH )
    continue;
IKReal x67=IKatan2(x66, x55);
j14array[0]=((x65)+(((-1.00000000000000)*(x67))));
sj14array[0]=IKsin(j14array[0]);
cj14array[0]=IKcos(j14array[0]);
j14array[1]=((3.14159265358979)+(((-1.00000000000000)*(x67)))+(((-1.00000000000000)*(x65))));
sj14array[1]=IKsin(j14array[1]);
cj14array[1]=IKcos(j14array[1]);
if( j14array[0] > IKPI )
{
    j14array[0]-=IK2PI;
}
else if( j14array[0] < -IKPI )
{    j14array[0]+=IK2PI;
}
j14valid[0] = true;
if( j14array[1] > IKPI )
{
    j14array[1]-=IK2PI;
}
else if( j14array[1] < -IKPI )
{    j14array[1]+=IK2PI;
}
j14valid[1] = true;
if( j14valid[0] && j14valid[1] && IKabs(cj14array[0]-cj14array[1]) < 0.0001 && IKabs(sj14array[0]-sj14array[1]) < 0.0001 )
{
    j14valid[1]=false;
}
for(int ij14 = 0; ij14 < 2; ++ij14)
{
if( !j14valid[ij14] )
{
    continue;
}
j14 = j14array[ij14]; cj14 = cj14array[ij14]; sj14 = sj14array[ij14];

IKReal soleval[1];
soleval[0]=((-0.0232000000000000)+(((-1.00000000000000)*(cj14)*(((0.0680000000000000)+(((-1.00000000000000)*(cj13)*(px)))+(((-1.00000000000000)*(py)*(sj13)))))))+(((-1.00000000000000)*(pz)*(sj14))));
if( soleval[0] > 0.0000000000000000  )
{
vsolutions.push_back(IKSolution()); IKSolution& solution = vsolutions.back();
solution.basesol.resize(2);
solution.basesol[0].foffset = j13;
solution.basesol[1].foffset = j14;
solution.vfree.resize(0);
}
}
}

}

}
}
}

}

}
}
return vsolutions.size()>0;
}

/// Durand-Kerner polynomial root finding method
static inline void polyroots4(IKReal rawcoeffs[4+1], IKReal rawroots[4], int& numroots)
{
    using std::complex;
    IKFAST_ASSERT(rawcoeffs[0] != 0);
    const IKReal tol = 128.0*std::numeric_limits<IKReal>::epsilon();
    complex<IKReal> coeffs[4];
    const int maxsteps = 50;
    for(int i = 0; i < 4; ++i) {
        coeffs[i] = complex<IKReal>(rawcoeffs[i+1]/rawcoeffs[0]);
    }
    complex<IKReal> roots[4];
    IKReal err[4];
    roots[0] = complex<IKReal>(1,0);
    roots[1] = complex<IKReal>(0.4,0.9); // any complex number not a root of unity is works
    err[0] = 1.0;
    err[1] = 1.0;
    for(int i = 2; i < 4; ++i) {
        roots[i] = roots[i-1]*roots[1];
        err[i] = 1.0;
    }
    for(int step = 0; step < maxsteps; ++step) {
        bool changed = false;
        for(int i = 0; i < 4; ++i) {
            if ( err[i] >= tol ) {
                changed = true;
                // evaluate
                complex<IKReal> x = roots[i] + coeffs[0];
                for(int j = 1; j < 4; ++j) {
                    x = roots[i] * x + coeffs[j];
                }
                for(int j = 0; j < 4; ++j) {
                    if( i != j ) {
                        if( roots[i] != roots[j] ) {
                            x /= (roots[i] - roots[j]);
                        }
                    }
                }
                roots[i] -= x;
                err[i] = abs(x);
            }
        }
        if( !changed ) {
            break;
        }
    }
    numroots = 0;
    for(int i = 0; i < 4; ++i) {
        if( IKabs(imag(roots[i])) < std::numeric_limits<IKReal>::epsilon() ) {
            rawroots[numroots++] = real(roots[i]);
        }
    }
}
};


/// solves the inverse kinematics equations.
/// \param pfree is an array specifying the free joints of the chain.
IKFAST_API bool ik(const IKReal* eetrans, const IKReal* eerot, const IKReal* pfree, std::vector<IKSolution>& vsolutions) {
IKSolver solver;
return solver.ik(eetrans,eerot,pfree,vsolutions);
}

IKFAST_API const char* getKinematicsHash() { return "80f514166e15c34bd64294fc1fdd5ddd"; }

#ifdef IKFAST_NAMESPACE
} // end namespace
#endif

#ifndef IKFAST_NO_MAIN
#include <stdio.h>
#include <stdlib.h>
#ifdef IKFAST_NAMESPACE
using namespace IKFAST_NAMESPACE;
#endif
int main(int argc, char** argv)
{
    if( argc != 12+getNumFreeParameters()+1 ) {
        printf("\nUsage: ./ik r00 r01 r02 t0 r10 r11 r12 t1 r20 r21 r22 t2 free0 ...\n\n"
               "Returns the ik solutions given the transformation of the end effector specified by\n"
               "a 3x3 rotation R (rXX), and a 3x1 translation (tX).\n"
               "There are %d free parameters that have to be specified.\n\n",getNumFreeParameters());
        return 1;
    }

    std::vector<IKSolution> vsolutions;
    std::vector<IKReal> vfree(getNumFreeParameters());
    IKReal eerot[9],eetrans[3];
    eerot[0] = atof(argv[1]); eerot[1] = atof(argv[2]); eerot[2] = atof(argv[3]); eetrans[0] = atof(argv[4]);
    eerot[3] = atof(argv[5]); eerot[4] = atof(argv[6]); eerot[5] = atof(argv[7]); eetrans[1] = atof(argv[8]);
    eerot[6] = atof(argv[9]); eerot[7] = atof(argv[10]); eerot[8] = atof(argv[11]); eetrans[2] = atof(argv[12]);
    for(std::size_t i = 0; i < vfree.size(); ++i)
        vfree[i] = atof(argv[13+i]);
    bool bSuccess = ik(eetrans, eerot, vfree.size() > 0 ? &vfree[0] : NULL, vsolutions);

    if( !bSuccess ) {
        fprintf(stderr,"Failed to get ik solution\n");
        return -1;
    }

    printf("Found %d ik solutions:\n", (int)vsolutions.size());
    std::vector<IKReal> sol(getNumJoints());
    for(std::size_t i = 0; i < vsolutions.size(); ++i) {
        printf("sol%d (free=%d): ", (int)i, (int)vsolutions[i].GetFree().size());
        std::vector<IKReal> vsolfree(vsolutions[i].GetFree().size());
        vsolutions[i].GetSolution(&sol[0],vsolfree.size()>0?&vsolfree[0]:NULL);
        for( std::size_t j = 0; j < sol.size(); ++j)
            printf("%.15f, ", sol[j]);
        printf("\n");
    }
    return 0;
}

#endif

#include "ikbase.h"
namespace IKFAST_NAMESPACE {
#ifdef RAVE_REGISTER_BOOST
#include BOOST_TYPEOF_INCREMENT_REGISTRATION_GROUP()
BOOST_TYPEOF_REGISTER_TYPE(IKSolution)
#endif
IkSolverBasePtr CreateIkSolver(EnvironmentBasePtr penv, const std::vector<dReal>& vfreeinc) {
    std::vector<int> vfree(getNumFreeParameters());
    for(size_t i = 0; i < vfree.size(); ++i) {
        vfree[i] = getFreeParameters()[i];
    }
    return IkSolverBasePtr(new IkFastSolver<IKReal,IKSolution>(ik,vfree,vfreeinc,getNumJoints(),static_cast<IkParameterizationType>(getIKType()), boost::shared_ptr<void>(), getKinematicsHash(), penv));
}
} // end namespace