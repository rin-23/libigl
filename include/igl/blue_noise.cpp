// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2020 Rinat Abdrashitov <rinat@dgp.toronto.edu>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.

#include "blue_noise.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <cassert>
#include <stdint.h>
#include <vector>
#include <Eigen/Core>

//TODO: reimplement bridsons code below using Eigen
namespace bridson
{
    template<class T>
    inline T sqr(const T &x)
    { return x*x; }

    // transforms even the sequence 0,1,2,3,... into reasonably good random numbers 
    // challenge: improve on this in speed and "randomness"!
    inline unsigned int randhash(unsigned int seed)
    {
        unsigned int i=(seed^12345391u)*2654435769u;
        i^=(i<<6)^(i>>26);
        i*=2654435769u;
        i+=(i<<5)^(i>>12);
        return i;
    }

    // returns repeatable stateless pseudo-random number in [0,1]
    inline double randhashd(unsigned int seed)
    { return randhash(seed)/(double)UINT_MAX; }
    inline float randhashf(unsigned int seed)
    { return randhash(seed)/(float)UINT_MAX; }

    // returns repeatable stateless pseudo-random number in [a,b]
    inline double randhashd(unsigned int seed, double a, double b)
    { return (b-a)*randhash(seed)/(double)UINT_MAX + a; }
    inline float randhashf(unsigned int seed, float a, float b)
    { return ( (b-a)*randhash(seed)/(float)UINT_MAX + a); }

    template<class T>
    void erase_unordered(std::vector<T> &a, unsigned int index)
    {
        a[index]=a.back();
        a.pop_back();
    }

    template<class T>
    void erase_unordered_swap(std::vector<T> &a, unsigned int index)
    {
        std::swap(a[index], a.back());
        a.pop_back();
    }

    template<class T>
    void find_and_erase_unordered(std::vector<T> &a, const T &doomed_element)
    {
        for(unsigned int i=0; i<a.size(); ++i)
            if(a[i]==doomed_element){
                erase_unordered(a, i);
                return;
            }
    }

    template<unsigned int N, class T>
    struct Vec
    {
        T v[N];

        Vec<N,T>(void)
        {}

        Vec<N,T>(T value_for_all)
        { for(unsigned int i=0; i<N; ++i) v[i]=value_for_all; }

        template<class S>
        Vec<N,T>(const S *source)
        { for(unsigned int i=0; i<N; ++i) v[i]=(T)source[i]; }

        Vec<N,T>(T v0, T v1)
        {
            assert(N==2);
            v[0]=v0; v[1]=v1;
        }

        Vec<N,T>(T v0, T v1, T v2)
        {
            assert(N==3);
            v[0]=v0; v[1]=v1; v[2]=v2;
        }

        T &operator[](int index)
        {
            assert(0<=index && (unsigned int)index<N);
            return v[index];
        }

        const T &operator[](int index) const
        {
            assert(0<=index && (unsigned int)index<N);
            return v[index];
        }

        Vec<N,T> operator+=(const Vec<N,T> &w)
        {
            for(unsigned int i=0; i<N; ++i) v[i]+=w[i];
            return *this;
        }

        Vec<N,T> operator+(const Vec<N,T> &w) const
        {
            Vec<N,T> sum(*this);
            sum+=w;
            return sum;
        }

        Vec<N,T> operator-=(const Vec<N,T> &w)
        {
            for(unsigned int i=0; i<N; ++i) v[i]-=w[i];
            return *this;
        }

        Vec<N,T> operator-(void) const // unary minus
        {
            Vec<N,T> negative;
            for(unsigned int i=0; i<N; ++i) negative.v[i]=-v[i];
            return negative;
        }

        Vec<N,T> operator-(const Vec<N,T> &w) const // (binary) subtraction
        {
            Vec<N,T> diff(*this);
            diff-=w;
            return diff;
        }

        Vec<N,T> operator*=(T a)
        {
            for(unsigned int i=0; i<N; ++i) v[i]*=a;
            return *this;
        }

        Vec<N,T> operator*(T a) const
        {
            Vec<N,T> w(*this);
            w*=a;
            return w;
        }

        Vec<N,T> operator*(const Vec<N,T> &w) const
        {
            Vec<N,T> componentwise_product;
            for(unsigned int i=0; i<N; ++i) componentwise_product[i]=v[i]*w.v[i];
            return componentwise_product;
        }

        Vec<N,T> operator/=(T a)
        {
            for(unsigned int i=0; i<N; ++i) v[i]/=a;
            return *this;
        }

        Vec<N,T> operator/(T a) const
        {
            Vec<N,T> w(*this);
            w/=a;
            return w;
        }
    };

    template<unsigned int N, class T>
    T mag2(const Vec<N,T> &a)
    {
        T l=sqr(a.v[0]);
        for(unsigned int i=1; i<N; ++i) l+=sqr(a.v[i]);
        return l;
    }

    template<unsigned int N, class T> 
    inline T dist2(const Vec<N,T> &a, const Vec<N,T> &b)
    { 
        T d=sqr(a.v[0]-b.v[0]);
        for(unsigned int i=1; i<N; ++i) d+=sqr(a.v[i]-b.v[i]);
        return d;
    }

    template<unsigned int N, class T> 
    inline bool operator==(const Vec<N,T> &a, const Vec<N,T> &b)
    { 
        bool t = (a.v[0] == b.v[0]);
        unsigned int i=1;
        while(i<N && t) {
            t = t && (a.v[i]==b.v[i]); 
            ++i;
        }
        return t;
    }

    template<unsigned int N, class T> 
    inline bool operator!=(const Vec<N,T> &a, const Vec<N,T> &b)
    { 
        bool t = (a.v[0] != b.v[0]);
        unsigned int i=1;
        while(i<N && !t) {
            t = t || (a.v[i]!=b.v[i]); 
            ++i;
        }
        return t;
    }

    template<unsigned int N, class T>
    inline Vec<N,T> operator*(T a, const Vec<N,T> &v)
    {
        Vec<N,T> w(v);
        w*=a;
        return w;
    }

    template<unsigned int N, class T>
    void sample_annulus(T radius, const Vec<N,T> &centre, unsigned int &seed, Vec<N,T> &x)
    {
        Vec<N,T> r;
        for(;;){
            for(unsigned int i=0; i<N; ++i){
                r[i]=4*(randhash(seed++)/(T)UINT_MAX-(T)0.5);
            }
            T r2=mag2(r);
            if(r2>1 && r2<=4)
                break;
        }
        x=centre+radius*r;
    }

    template<unsigned int N, class T>
    unsigned long int n_dimensional_array_index(const Vec<N,unsigned int> &dimensions, const Vec<N,T> &x)
    {
        unsigned long int k=0;
        if(x[N-1]>=0){
            k=(unsigned long int)x[N-1];
            if(k>=dimensions[N-1]) k=dimensions[N-1]-1;
        }
        for(unsigned int i=N-1; i>0; --i){
            k*=dimensions[i-1];
            if(x[i-1]>=0){
                unsigned int j=(int)x[i-1];
                if(j>=dimensions[i-1]) j=dimensions[i-1]-1;
                k+=j;
            }
        }
        return k;
    }

    template<unsigned int N, class T>
    void bluenoise_sample(T radius, Vec<N,T> xmin, Vec<N,T> xmax, std::vector<Vec<N,T> > &sample,
                        unsigned int seed=0, int max_sample_attempts=30)
    {
        sample.clear();
        std::vector<unsigned int> active_list;

        // acceleration grid
        T grid_dx=T(0.999)*radius/std::sqrt((T)N); // a grid cell this size can have at most one sample in it
        Vec<N,unsigned int> dimensions;
        unsigned long int total_array_size=1;
        for(unsigned int i=0; i<N; ++i){
            dimensions[i]=(unsigned int)std::ceil((xmax[i]-xmin[i])/grid_dx);
            total_array_size*=dimensions[i];
        }
        std::vector<int> accel(total_array_size, -1); // -1 indicates no sample there; otherwise index of sample point

        // first sample
        Vec<N,T> x;
        for(unsigned int i=0; i<N; ++i){
            x[i]=(xmax[i]-xmin[i])*(randhash(seed++)/(T)UINT_MAX) + xmin[i];
        }
        sample.push_back(x);
        active_list.push_back(0);
        unsigned int k=n_dimensional_array_index(dimensions, (x-xmin)/grid_dx);
        accel[k]=0;

        while(!active_list.empty()){
            unsigned int r=int(randhashf(seed++, 0, active_list.size()-0.0001f));
            int p=active_list[r];
            bool found_sample=false;
            Vec<N,unsigned int> j, jmin, jmax;
            for(int attempt=0; attempt<max_sample_attempts; ++attempt){
                sample_annulus(radius, sample[p], seed, x);
                // check this sample is within bounds
                for(unsigned int i=0; i<N; ++i){
                    if(x[i]<xmin[i] || x[i]>xmax[i])
                    goto reject_sample;
                }
                // test proximity to nearby samples
                for(unsigned int i=0; i<N; ++i){
                    int thismin=(int)((x[i]-radius-xmin[i])/grid_dx);
                    if(thismin<0) thismin=0;
                    else if(thismin>=(int)dimensions[i]) thismin=dimensions[i]-1;
                    jmin[i]=(unsigned int)thismin;
                    int thismax=(int)((x[i]+radius-xmin[i])/grid_dx);
                    if(thismax<0) thismax=0;
                    else if(thismax>=(int)dimensions[i]) thismax=dimensions[i]-1;
                    jmax[i]=(unsigned int)thismax;
                }
                for(j=jmin;;){
                    // check if there's a sample at j that's too close to x
                    k=n_dimensional_array_index(dimensions, j);
                    if(accel[k]>=0 && accel[k]!=p){ // if there is a sample point different from p
                        if(dist2(x, sample[accel[k]])<sqr(radius))
                            goto reject_sample;
                        }
                        // move on to next j
                        for(unsigned int i=0; i<N; ++i){
                        ++j[i];
                        if(j[i]<=jmax[i]){
                            break;
                        }else{
                            if(i==N-1) goto done_j_loop;
                            else j[i]=jmin[i]; // and try incrementing the next dimension along
                        }
                    }
                }
                done_j_loop:
                // if we made it here, we're good!
                found_sample=true;
                break;
                // if we goto here, x is too close to an existing sample
                reject_sample:
                ; // nothing to do except go to the next iteration in this loop
            }
            if(found_sample){
                size_t q=sample.size(); // the index of the new sample
                sample.push_back(x);
                active_list.push_back(q);
                k=n_dimensional_array_index(dimensions, (x-xmin)/grid_dx);
                accel[k]=(int)q;
            }else{
                // since we couldn't find a sample on p's disk, we remove p from the active list
                erase_unordered(active_list, r);
            }
        }
    }
}

template <unsigned int N, typename DerivedX, typename DerivedS>
IGL_INLINE void igl::blue_noise(typename DerivedS::Scalar radius, 
                                const Eigen::PlainObjectBase<DerivedX>& xmin,  
                                const Eigen::PlainObjectBase<DerivedX>& xmax, 
                                unsigned int seed, 
                                int max_sample_attempts, 
                                Eigen::PlainObjectBase<DerivedS>& S)
{
    assert(N > 0 && "dimension must be positive");
    assert(radius > 0 && "radius must be positive");
    assert(xmin.size() == N &&  xmax.size() == N && "invalid extent of the sample domain");

    if (N <= 0 || radius <= 0 || xmin.size() != N || xmax.size() != N)
        return;

    typedef typename DerivedS::Scalar Scalar;
    using namespace bridson;

    Vec<N,Scalar> xmin2;
    Vec<N,Scalar> xmax2;
    for (int i=0; i < N; ++i) 
    {
        xmin2[i] = xmin(i);
        xmax2[i] = xmax(i);
    }

    std::vector<Vec<N,Scalar>> sample;
    bluenoise_sample(radius, xmin2, xmax2, sample, seed, max_sample_attempts);

    S.resize(sample.size(), N);
    for (int i = 0; i < sample.size(); ++i) 
        for (int j = 0; j < N; ++j) 
            S(i,j) = sample[i][j];
}
