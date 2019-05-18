//============================================================================
// Name        : Polynomial.cpp
// Author      : Rebecca Su
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <cmath>
#include "Polynomial.h"

#ifndef MARMOSET_TESTING
int main();
#endif

void point_to_new(poly_t &p, double * p_new);

#ifndef MARMOSET_TESTING
int main()
{
	unsigned int test_degree{1};
	double a_test[2]{1, 1};
	poly_t p_test;
	poly_t q_test;

	init_poly(p_test, a_test, test_degree);
	init_poly(q_test, a_test, test_degree);

	poly_add(p_test, q_test);
	std::cout << p_test.degree;
	destroy_poly(p_test);
	destroy_poly(q_test);

	return 0;
}
#endif

void point_to_new(poly_t &p, double * p_new)
{
	delete[] p.a_coeffs;
	p.a_coeffs = p_new;
}

void init_poly(poly_t &p, double const init_coeffs[], unsigned int const init_degree)
{
	if(p.a_coeffs != nullptr)
	{
		delete[] p.a_coeffs;
		p.a_coeffs = nullptr;
	}

	p.a_coeffs = new double[init_degree + 1];
	p.degree = init_degree;
	for(int k = 0; k <= init_degree; k++)
		p.a_coeffs[k] = init_coeffs[k];
}

void destroy_poly(poly_t &p)
{
	delete[] p.a_coeffs;
	p.a_coeffs = nullptr;
}

unsigned int poly_degree(poly_t const &p)
{
	if(p.a_coeffs == nullptr)
		throw 0;
	else
		return p.degree;
}

double poly_coeff(poly_t const &p, unsigned int n)
{
	double coefficient{0};

	if(p.a_coeffs == nullptr)
		throw 0;
	else
	{
		if(n <= p.degree)
			coefficient = p.a_coeffs[n];
		else
			coefficient = 0;
	}
	return coefficient;
}

double poly_val(poly_t const &p, double const x)
{
	double value{0};

	if(p.a_coeffs == nullptr)
		throw 0;
	else
	{
		for(int k = 0; k < p.degree + 1; k ++)
			value += pow(x, k) * p.a_coeffs[k];
	}
	return value;
}

void poly_add(poly_t &p, poly_t const &q)
{
	double *p_newp;
	bool is_zero{true};

	if(p.a_coeffs == nullptr || q.a_coeffs == nullptr)
		throw 0;
	else
	{
		p_newp = new double [std::max(p.degree, q.degree) + 1];
		for(int k = 0; k < std::min(p.degree, q.degree) + 1; k++)
			p_newp[k] = p.a_coeffs[k] + q.a_coeffs[k];
		for(int k = std::min(p.degree, q.degree) + 1; k < std::max(p.degree, q.degree) + 1; k++)
		{
			if(p.degree >= q.degree)
				p_newp[k] = p.a_coeffs[k];
			else
				p_newp[k] = q.a_coeffs[k];
		}
	}

	point_to_new(p, p_newp);

	p.degree = std::max(p.degree, q.degree);

	while(is_zero == true)
	{
		for(int k = p.degree; k >= 0; k--)
			{
				if(p.a_coeffs[k] == 0)
					p.degree--;
				else
					is_zero = false;
			}
	}

}

void poly_subtract(poly_t &p, poly_t const &q)
{
	double *p_newp;
	bool is_zero{true};

	if(p.a_coeffs == nullptr || q.a_coeffs == nullptr)
		throw 0;
	else
	{
		p_newp = new double [std::max(p.degree, q.degree) + 1];
		for(int k = 0; k < std::min(p.degree, q.degree) + 1; k++)
			p_newp[k] = p.a_coeffs[k] - q.a_coeffs[k];
		for(int k = std::min(p.degree, q.degree) + 1; k < std::max(p.degree, q.degree) + 1; k++)
		{
			if(p.degree >= q.degree)
				p_newp[k] = p.a_coeffs[k];
			else
				p_newp[k] = -q.a_coeffs[k];
		}
	}

	point_to_new(p, p_newp);

	p.degree = std::max(p.degree, q.degree);

	while(is_zero == true)
	{
		for(int k = p.degree; k >= 0; k--)
			{
				if(p.a_coeffs[k] == 0)
					p.degree--;
				else
					is_zero = false;
			}
	}

}

void poly_multiply(poly_t &p, poly_t const &q)
{
	double * p_new;
	unsigned int k{p.degree + q.degree};

	if(p.a_coeffs == nullptr || q.a_coeffs == nullptr)
		throw 0;
	else
	{
		p_new = new double [k + 2];
		for(int k = 0; k <= p.degree + q.degree; k++)
			p_new[k] = 0;

		for(int x = 0; x < p.degree + 1; x++)
		{
			for(int y = 0; y < q.degree + 1; y++)
			{
				p_new[x + y] += p.a_coeffs[x] * q.a_coeffs[y];
			}
		}
	}
	point_to_new(p, p_new);
	p.degree = p.degree + q.degree;
}

double poly_divide(poly_t &p, double r)
{
	double remainder{0};

	if(p.a_coeffs == nullptr)
		throw 0;

	if(p.degree == 0)
	{
		remainder = p.a_coeffs[0];
		p.a_coeffs[0] = 0;
	}
	else
	{
		double * p_new = new double[p.degree + 1];
		remainder = p.a_coeffs[p.degree];

		for(int k = p.degree; k >= 1; k--)
		{
				p_new[k - 1] = remainder;
				remainder = r * remainder + p.a_coeffs[k - 1];
		}
		point_to_new(p, p_new);
		p.degree = p.degree - 1;
	}
	return remainder;
}

void poly_diff(poly_t &p)
{
	double *p_new{nullptr};
	if(p.a_coeffs == nullptr)
		throw 0;
	else if(p.degree == 0)
	{
		p.a_coeffs = new double[1]{0};
		p.degree = 0;
	}else
	{
		p_new = new double[p.degree + 1];
		for(int k = 0; k <= p.degree; k++)
			p_new[k] = p.a_coeffs[k + 1] * (k + 1);
		p.degree = p.degree - 1;
	}
	point_to_new(p, p_new);
}

double poly_approx_int(poly_t const &p, double a, double b, unsigned int n)
{
	double h{0};
	double sum{0};

	if(p.a_coeffs == nullptr)
		throw 0;
	else
	{
		h = (b - a) / n;
		for(int k = 1; k < n; k++)
			sum += 2 * poly_val(p, a + k * h);
		sum = h / 2 * (sum + poly_val(p, a) + poly_val(p, b));
		return sum;
	}
}
