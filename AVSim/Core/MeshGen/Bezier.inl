#ifndef MESHGEN_BEZIER_INL
#define MESHGEN_BEZIER_INL

namespace World3d{
    template<int N, typename Field>
    DReal ArcLength(const BezierCurve<N, Vector>& x, DReal t_st, DReal t_end, DReal err){
        QuickLineIntegral qli;
        qli.m_curve = [&x](double t) { return CGAL::ORIGIN + x.eval(t); };
        qli.Tst=t_st, qli.Tend=t_end; qli.Ninit = (N + 1)*10; qli.rel=err;
        qli.setup();
        return qli(t_end);
    }

    template<int N, int K, typename Field, typename MakeFieldProduct>
    BezierCurve<N+K, decltype(std::declval<MakeFieldProduct>().operator()(std::declval<Field>(), std::declval<Field>()))>
        PerformBezierProduct(const BezierCurve<N, Field>& a, const BezierCurve<K, Field>& b, MakeFieldProduct prod){
        using ProdType = decltype(std::declval<MakeFieldProduct>().operator()(std::declval<Field>(), std::declval<Field>()));
        BezierCurve<N+K, ProdType> r;
        std::fill(r.B+1, r.B + N + K, ProdType());
        r.B[0] = a.B[0] * b.B[0];
        r.B[N + K] = a.B[N] * b.B[K];
        long long int c_k = 1;
        for (int k = 1; k < N + K; ++k){
            c_k *= (N + K + 1 - k);
            c_k /= k;
            int I_min = std::max(0, k-K), I_max = std::min(k, N);
            long long int c_i = 1, c_j = 1;
            for (int i = 1; i < I_min; ++i){
                c_i *= (N + 1 - i);
                c_i /= i;
            }
            int j_st = std::min(K-k + I_min-1, k - I_min + 1);
            for (int i = 1; i < j_st; ++i){
                c_i *= (K + 1 - i);
                c_i /= i;
            }
            for (int i = I_min; i < I_max + 1; ++i){
                int j = k - i;
                if (i > 0){
                    c_i *= (N + 1 - i); 
                    c_i /= i;
                }
                if (K - j > 0){
                    c_j *= (j+1);
                    c_j /= (K - j);
                }
                r.B[k] += c_i * c_j * prod(a.B[i], b.B[j]) / c_k;  
            }
        }
        return r;
    }

    template<int N, int K, typename Field>
    BezierCurve<N+K, decltype(std::declval<Field>()*std::declval<Field>())> operator*(const BezierCurve<N, Field>& a, const BezierCurve<K, Field>& b){
        using ProdType = decltype(std::declval<Field>()*std::declval<Field>());
        auto prod = [](const Field& a, const Field& b) -> ProdType { return a*b; };
        return PerformBezierProduct<N, K, Field>(a, b, prod);
    }

    template<int N, typename Field> template<int K, typename std::enable_if< (K>N) >::type*>
    BezierCurve<K, Field> BezierCurve<N, Field>::raise_order() const{
        BezierCurve<K, Field> r;
        std::fill(r.B+1, r.B + K, Field());
        r.B[0] = B[0];
        r.B[K] = B[N];
        long long int c_k = 1;
        for (int k = 1; k < K; ++k){
            c_k *= (K + 1 - k);
            c_k /= k; 
            int I_min = std::max(0, k-(K-N)), I_max = std::min(k, N);
            long long int c_i = 1, c_j = 1;
            for (int i = 1; i < I_min; ++i){
                c_i *= (N + 1 - i);
                c_i /= i;
            }
            int j_st = std::min((K-N)-k + I_min-1, k - I_min + 1);
            for (int i = 1; i < j_st; ++i){
                c_i *= ((K-N) + 1 - i);
                c_i /= i;
            }
            for (int i = I_min; i < I_max + 1; ++i){
                int j = k - i;
                if (i > 0){
                    c_i *= (N + 1 - i); 
                    c_i /= i;
                }
                if ((K-N) - j > 0){
                    c_j *= (j+1);
                    c_j /= ((K-N) - j);
                }
                r.B[k] += c_i * c_j * B[i] / c_k;  
            }
        }
        return r;
    }

    template<int N, typename Field>
    std::array<Field, N+1> BezierCurve<N, Field>::polynomial() const{
        std::array<Field, N+1> p{0};
        long long int c_j = 1;
        for (int j = 0; j < N+1; ++j){
            if (j > 0) {
                c_j *= (N+1-j);
                c_j /= j;
            }
            long long int c_i = 1;
            for (int i = 0; i < j+1; ++i){
                if (i > 0){
                    c_i *= j + 1 - i;
                    c_i /= i;
                }
                if (((i + j) % 2) == 0)
                    p[j] += c_i * c_j * B[i];
                else 
                    p[j] -= c_i * c_j * B[i];
            }
        }
        return p;
    }
    template<int N, typename Field>
    BezierCurve<N-1, Field> BezierCurve<N, Field>::derivative() const{
        BezierCurve<N-1, Field> b;
        for (int i = 0; i < N; ++i) b.B[i] = N*(B[i+1]-B[i]);
        return b;
    }
    template<int N, typename Field>
    Field BezierCurve<N, Field>::eval(DReal t) const{
        Field tmp[N];
        auto q = (1 - t);
        for (int k = 0; k < N; ++k) 
            tmp[k] = B[k] * q + B[k + 1] * t;
        
        for (int j = 1; j < N; ++j) {
            for (int k = 0; k < N - j; ++k) {
                tmp[k] = tmp[k] * q + tmp[k + 1] * t;
            }
        }
        return tmp[0];
    }
    template<int N, typename Field>
    std::pair<BezierCurve<N, Field>, BezierCurve<N, Field>> BezierCurve<N, Field>::splitAt(DReal t_split) const{
        BezierCurve<N, Field> b1, b2;
        b1.B[0] = B[0];
        std::copy(B, B + N + 1, b2.B);
        auto q = (1 - t_split);


        for (int j = 1; j < N+1; ++j) {
            for (int k = 0; k < N+1 - j; ++k) {
                b2.B[k] = b2.B[k] * q + b2.B[k + 1] * t_split;
            }
            b1.B[j] = b2.B[0];
        }

        return {std::move(b1), std::move(b2)};
    }

    template<int N, typename Field1, typename Field2> 
        BezierCurve<N, decltype(declval<Field1>()*declval<Field2>())> operator*(const Field1& a, const BezierCurve<N, Field2>& b){
        using ProdType = decltype(declval<Field1>()*declval<Field2>());
        std::array<ProdType, N+1> B = {a*b.B[0]};
        for (int i = 0; i < N+1; ++i) B[i] = a*b.B[i];
        return {B};
    }
    template<int N, typename Field1, typename Field2> 
        BezierCurve<N, decltype(declval<Field1>()*declval<Field2>())> operator*(const BezierCurve<N, Field1>& a, const Field2& b){
        using ProdType = decltype(declval<Field1>()*declval<Field2>());
        std::array<ProdType, N+1> B = {a.B[0]*b};
        for (int i = 0; i < N+1; ++i) B[i] = a.B[i]*b;
        return {B};
    }
    template<int N, typename Field> std::ostream& operator<<(std::ostream& out, const BezierCurve<N, Field>& a){
        out << "BezierCurve<" << N << ">:{ ";
        for (int i = 0; i < N+1; ++i) out << "(" << a.B[i] << ") ";
        out << "}";
        return out;
    }
}

#endif