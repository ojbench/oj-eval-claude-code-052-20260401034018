#include <iostream>
#include <algorithm>
#include <cstring>
#include <string>
#include <cctype>

class term {
public:
    int a, b, c, d;

    term() : a(0), b(0), c(0), d(0) {}
    term(int _a, int _b, int _c, int _d) : a(_a), b(_b), c(_c), d(_d) {}

    bool operator == (const term &obj) const {
        return b == obj.b && c == obj.c && d == obj.d;
    }
    bool operator != (const term &obj) const {
        return b != obj.b || c != obj.c || d != obj.d;
    }
    bool operator < (const term &obj) const {
        if (b != obj.b) return b > obj.b;
        if (c != obj.c) return c > obj.c;
        return d > obj.d;
    }
};

class poly {
public:
    int n;
    term *t;

    poly() : n(0), t(NULL) {}
    poly(int _n) {
        n = _n;
        t = new term[n];
    }
    poly(const poly &p) {
        n = p.n;
        t = new term[n];
        for (int i = 0; i < n; ++i) {
            t[i] = p.t[i];
        }
    }

    void simplify() {
        if (n == 0) return;
        std::sort(t, t + n);

        int writeIdx = 0;
        for (int i = 0; i < n; ++i) {
            if (t[i].a == 0) continue;

            if (writeIdx > 0 && t[writeIdx-1] == t[i]) {
                t[writeIdx-1].a += t[i].a;
                if (t[writeIdx-1].a == 0) writeIdx--;
            } else {
                t[writeIdx++] = t[i];
            }
        }
        n = writeIdx;
    }

    poly operator + (const poly &obj) const {
        poly ans(n + obj.n);
        for (int i = 0; i < n; ++i) ans.t[i] = t[i];
        for (int i = 0; i < obj.n; ++i) ans.t[i + n] = obj.t[i];
        ans.simplify();
        return ans;
    }

    poly operator - (const poly &obj) const {
        poly ans(n + obj.n);
        for (int i = 0; i < n; ++i) ans.t[i] = t[i];
        for (int i = 0; i < obj.n; ++i) {
            ans.t[i + n] = obj.t[i];
            ans.t[i + n].a *= -1;
        }
        ans.simplify();
        return ans;
    }

    poly operator * (const poly &obj) const {
        poly ans(n * obj.n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < obj.n; ++j) {
                int idx = i * obj.n + j;
                ans.t[idx].a = t[i].a * obj.t[j].a;
                ans.t[idx].b = t[i].b + obj.t[j].b;
                ans.t[idx].c = t[i].c + obj.t[j].c;
                ans.t[idx].d = t[i].d + obj.t[j].d;
            }
        }
        ans.simplify();
        return ans;
    }

    poly& operator = (const poly &obj) {
        if (&obj == this) return *this;
        n = obj.n;
        if (t != NULL) delete []t;
        t = new term[n];
        for (int i = 0; i < n; ++i) {
            t[i] = obj.t[i];
        }
        return *this;
    }

    poly derivate() const {
        poly ans(n * 3); // At most 3 terms per term
        int idx = 0;
        for (int i = 0; i < n; ++i) {
            // Derivative of x^b term
            if (t[i].b > 0) {
                ans.t[idx] = term(t[i].a * t[i].b, t[i].b - 1, t[i].c, t[i].d);
                idx++;
            }
            // Derivative of sin^c x term: c*sin^(c-1)x * cosx
            if (t[i].c > 0) {
                ans.t[idx] = term(t[i].a * t[i].c, t[i].b, t[i].c - 1, t[i].d + 1);
                idx++;
            }
            // Derivative of cos^d x term: -d*cos^(d-1)x * sinx
            if (t[i].d > 0) {
                ans.t[idx] = term(-t[i].a * t[i].d, t[i].b, t[i].c + 1, t[i].d - 1);
                idx++;
            }
        }
        ans.n = idx;
        ans.simplify();
        return ans;
    }

    ~poly() {
        if (t != NULL) delete []t;
    }
};

class frac {
public:
    poly p, q;

    frac() {}
    frac(int x) {
        p = poly(1);
        p.t[0] = term(x, 0, 0, 0);
        q = poly(1);
        q.t[0] = term(1, 0, 0, 0);
    }
    frac(term _p) {
        q = poly(1);
        q.t[0] = term(1, 0, 0, 0);
        p = poly(1);
        p.t[0] = _p;
    }
    frac(poly _p, poly _q) : p(_p), q(_q) {}

    frac operator + (const frac &obj) const {
        return frac(p * obj.q + q * obj.p, q * obj.q);
    }
    frac operator - (const frac &obj) const {
        return frac(p * obj.q - q * obj.p, q * obj.q);
    }
    frac operator * (const frac &obj) const {
        return frac(p * obj.p, q * obj.q);
    }
    frac operator / (const frac &obj) const {
        return frac(p * obj.q, q * obj.p);
    }

    frac derivate() const {
        poly dp = p.derivate();
        poly dq = q.derivate();
        return frac(dp * q - dq * p, q * q);
    }

    void output() {
        // Check if numerator is zero
        if (p.n == 0) {
            std::cout << "0" << std::endl;
            return;
        }

        // Check if denominator is 1
        bool denomIsOne = (q.n == 1 && q.t[0].a == 1 && q.t[0].b == 0 &&
                          q.t[0].c == 0 && q.t[0].d == 0);

        if (denomIsOne) {
            outputPoly(p, false);
            std::cout << std::endl;
        } else {
            // Output as (numerator)/(denominator)
            bool needParenP = (p.n > 1);
            bool needParenQ = (q.n > 1);

            if (needParenP) std::cout << "(";
            outputPoly(p, false);
            if (needParenP) std::cout << ")";

            std::cout << "/";

            if (needParenQ) std::cout << "(";
            outputPoly(q, false);
            if (needParenQ) std::cout << ")";

            std::cout << std::endl;
        }
    }

    void outputPoly(const poly &pol, bool) {
        for (int i = 0; i < pol.n; ++i) {
            const term &tm = pol.t[i];
            int coef = tm.a;

            if (i == 0) {
                if (coef < 0) {
                    std::cout << "-";
                    coef = -coef;
                }
            } else {
                if (coef < 0) {
                    std::cout << "-";
                    coef = -coef;
                } else {
                    std::cout << "+";
                }
            }

            bool isConstant = (tm.b == 0 && tm.c == 0 && tm.d == 0);

            if (coef != 1 || isConstant) {
                std::cout << coef;
            }

            if (tm.b > 0) {
                std::cout << "x";
                if (tm.b > 1) std::cout << "^" << tm.b;
            }

            if (tm.c > 0) {
                std::cout << "sin";
                if (tm.c > 1) std::cout << "^" << tm.c;
                std::cout << "x";
            }

            if (tm.d > 0) {
                std::cout << "cos";
                if (tm.d > 1) std::cout << "^" << tm.d;
                std::cout << "x";
            }
        }
    }
};

// Parser
int pos;
frac parseExpr(const char *s);
frac parseTerm(const char *s);
frac parseAtom(const char *s);

frac parseAtom(const char *s) {
    if (s[pos] == '(') {
        pos++;
        frac result = parseExpr(s);
        pos++; // skip ')'
        return result;
    }

    // Parse a basic term: coefficient * x^exp * sin^exp * cos^exp
    int sign = 1;
    if (s[pos] == '-') {
        sign = -1;
        pos++;
    } else if (s[pos] == '+') {
        pos++;
    }

    int coef = 0;
    int xExp = 0, sinExp = 0, cosExp = 0;

    // Parse coefficient (if any)
    if (isdigit(s[pos])) {
        while (isdigit(s[pos])) {
            coef = coef * 10 + (s[pos] - '0');
            pos++;
        }
    } else {
        coef = 1;
    }

    // Parse x, sinx, cosx
    while (s[pos] && s[pos] != '+' && s[pos] != '-' && s[pos] != '*' &&
           s[pos] != '/' && s[pos] != ')') {
        if (s[pos] == 'x') {
            pos++;
            if (s[pos] == '^') {
                pos++;
                int exp = 0;
                while (isdigit(s[pos])) {
                    exp = exp * 10 + (s[pos] - '0');
                    pos++;
                }
                xExp = exp;
            } else {
                xExp = 1;
            }
        } else if (s[pos] == 's') { // sinx
            pos += 3; // skip "sin"
            if (s[pos] == '^') {
                pos++;
                int exp = 0;
                while (isdigit(s[pos])) {
                    exp = exp * 10 + (s[pos] - '0');
                    pos++;
                }
                pos++; // skip 'x'
                sinExp = exp;
            } else {
                pos++; // skip 'x'
                sinExp = 1;
            }
        } else if (s[pos] == 'c') { // cosx
            pos += 3; // skip "cos"
            if (s[pos] == '^') {
                pos++;
                int exp = 0;
                while (isdigit(s[pos])) {
                    exp = exp * 10 + (s[pos] - '0');
                    pos++;
                }
                pos++; // skip 'x'
                cosExp = exp;
            } else {
                pos++; // skip 'x'
                cosExp = 1;
            }
        } else {
            break;
        }
    }

    return frac(term(sign * coef, xExp, sinExp, cosExp));
}

frac parseTerm(const char *s) {
    frac result = parseAtom(s);

    while (s[pos] == '*' || s[pos] == '/') {
        char op = s[pos];
        pos++;
        frac right = parseAtom(s);
        if (op == '*') {
            result = result * right;
        } else {
            result = result / right;
        }
    }

    return result;
}

frac parseExpr(const char *s) {
    frac result = parseTerm(s);

    while (s[pos] == '+' || (s[pos] == '-' && pos > 0)) {
        char op = s[pos];
        pos++;
        frac right = parseTerm(s);
        if (op == '+') {
            result = result + right;
        } else {
            result = result - right;
        }
    }

    return result;
}

void solve(char *s, int n) {
    pos = 0;
    frac result = parseExpr(s);
    result.output();
    frac derivative = result.derivate();
    derivative.output();
}

int main() {
    std::string str;
    std::cin >> str;
    int n = str.length();
    char *s = new char[n + 2]{0};
    for (int i = 0; i < n; ++i) s[i] = str[i];
    solve(s, n);
    delete []s;
    return 0;
}
