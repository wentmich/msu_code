# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

//# include "filon.h"

/******************************************************************************/

double filon_fun_cos ( int n, double *f ( int n, double x[] ), double a, 
  double b, double t )

/******************************************************************************/
/*
  Purpose:

    FILON_FUN_COS uses Filon's method on integrals with a cosine factor.

  Discussion:

    The integral to be approximated has the form:

      Integral ( A <= X <= B ) F(X) * COS(T*X) dX

    where T is user specified.

    The function is interpolated over each subinterval by
    a parabolic arc.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 May 2014

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

    Stephen Chase, Lloyd Fosdick,
    An Algorithm for Filon Quadrature,
    Communications of the Association for Computing Machinery,
    Volume 12, Number 8, August 1969, pages 453-457.

    Stephen Chase, Lloyd Fosdick,
    Algorithm 353:
    Filon Quadrature,
    Communications of the Association for Computing Machinery,
    Volume 12, Number 8, August 1969, pages 457-458.

    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.

  Parameters:

    Input, int N, the number of data points.
    N must be odd, and greater than 1.

    Input, double *F ( int n, double x[] ), the function which evaluates the 
    integrand.

    Input, double A, B, the limits of integration.

    Input, double T, the multiplier of the X argument of the cosine.

    Output, double FILON_FUN_COS, the approximate value of the integral.
*/
{
  double alpha;
  double beta;
  double c2n;
  double c2nm1;
  double cost;
  double *ftab;
  double gamma;
  double h;
  int i;
  double sint;
  double theta;
  double value;
  double *x;

  if ( a == b )
  {
    value = 0.0;
    return value;
  }
 
  if ( n <= 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "FILON_FUN_COS - Fatal error!\n" );
    fprintf ( stderr, "  N < 2\n" );
    fprintf ( stderr, "  N = %d\n", n );
    exit ( 1 );
  }
 
  if ( ( n % 2 ) != 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "FILON_FUN_COS - Fatal error!\n" );
    fprintf ( stderr, "  N must be odd.\n" );
    fprintf ( stderr, "  N = %d\n", n );
    exit ( 1 );
  }
/*
  Set the X values.
*/
  x = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( n - i - 1 ) * a   
           + ( double ) (     i     ) * b ) 
           / ( double ) ( n     - 1 );
  }

  h = ( b - a ) / ( double ) ( n - 1 );
  theta = t * h;
  sint = sin ( theta );
  cost = cos ( theta );

  if ( 6.0 * fabs ( theta ) <= 1.0 )
  {
    alpha = 2.0 * pow ( theta, 3 ) /   45.0 
          - 2.0 * pow ( theta, 5 ) /  315.0 
          + 2.0 * pow ( theta, 7 ) / 4725.0;
  
    beta =  2.0                    /     3.0 
          + 2.0 * pow ( theta, 2 ) /    15.0 
          - 4.0 * pow ( theta, 4 ) /   105.0 
          + 2.0 * pow ( theta, 6 ) /   567.0 
          - 4.0 * pow ( theta, 8 ) / 22275.0;

    gamma = 4.0                    /      3.0 
          - 2.0 * pow ( theta, 2 ) /     15.0 
          +       pow ( theta, 4 ) /    210.0 
          -       pow ( theta, 6 ) /  11340.0;
  }
  else
  {
    alpha = ( pow ( theta, 2 ) + theta * sint * cost - 2.0 * sint * sint ) 
      / pow ( theta, 3 );

    beta = ( 2.0 * theta + 2.0 * theta * cost * cost
      - 4.0 * sint * cost ) / pow ( theta, 3 );

    gamma = 4.0 * ( sint - theta * cost ) / pow ( theta, 3 );
  }
/*
  Tabulate the function.
*/
  ftab = f ( n, x );

  c2n = 0.5 * ftab[0] * cos ( t * x[0] );
  for ( i = 2; i < n - 1; i = i + 2 )
  {
    c2n = c2n + ftab[i] * cos ( t * x[i] );
  }
  c2n = c2n + 0.5 * ftab[n-1] * cos ( t * x[n-1] );

  c2nm1 = 0.0;
  for ( i = 1; i <= n - 2; i = i + 2 )
  {
    c2nm1 = c2nm1 + ftab[i] * cos ( t * x[i] );
  }

  value = h * ( 
      alpha * ( ftab[n-1] * sin ( t * x[n-1] )  
              - ftab[0]   * sin ( t * x[0] ) ) 
    + beta * c2n 
    + gamma * c2nm1 );

  free ( ftab );
  free ( x );

  return value;
}
/******************************************************************************/

double filon_tab_cos ( int n, double ftab[], double a, double b, double t )

/******************************************************************************/
/*
  Purpose:

    FILON_TAB_COS uses Filon's method on integrals with a cosine factor.

  Discussion:

    The integral to be approximated has the form:

      Integral ( A <= X <= B ) F(X) * COS(T*X) dX

    where T is user specified.

    The function is interpolated over each subinterval by
    a parabolic arc.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 May 2014

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

    Stephen Chase, Lloyd Fosdick,
    An Algorithm for Filon Quadrature,
    Communications of the Association for Computing Machinery,
    Volume 12, Number 8, August 1969, pages 453-457.

    Stephen Chase, Lloyd Fosdick,
    Algorithm 353:
    Filon Quadrature,
    Communications of the Association for Computing Machinery,
    Volume 12, Number 8, August 1969, pages 457-458.

    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.

  Parameters:

    Input, int N, the number of data points.
    N must be odd, and greater than 1.

    Input, double FTAB[N], contains the value of the function
    at A, A+H, A+2*H, ... , B-H, B, where H = (B-A)/(N-1).

    Input, double A, B, the limits of integration.

    Input, double T, the multiplier of the X argument of the cosine.

    Output, double FILON_TAB_COS, the approximate value of the integral.
*/
{
  double alpha;
  double beta;
  double c2n;
  double c2nm1;
  double cost;
  double gamma;
  double h;
  int i;
  double sint;
  double theta;
  double value;
  double *x;

  if ( a == b )
  {
    value = 0.0;
    return value;
  }
 
  if ( n <= 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "FILON_TAB_COS - Fatal error!\n" );
    fprintf ( stderr, "  N < 2\n" );
    fprintf ( stderr, "  N = %d\n", n );
    exit ( 1 );
  }
 
  if ( ( n % 2 ) != 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "FILON_TAB_COS - Fatal error!\n" );
    fprintf ( stderr, "  N must be odd.\n" );
    fprintf ( stderr, "  N = %d\n", n );
    exit ( 1 );
  }
/*
  Set the X values.
*/
  x = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( n - i - 1 ) * a   
           + ( double ) (     i     ) * b ) 
           / ( double ) ( n     - 1 );
  }

  h = ( b - a ) / ( double ) ( n - 1 );
  theta = t * h;
  sint = sin ( theta );
  cost = cos ( theta );

  if ( 6.0 * fabs ( theta ) <= 1.0 )
  {
    alpha = 2.0 * pow ( theta, 3 ) /   45.0 
          - 2.0 * pow ( theta, 5 ) /  315.0 
          + 2.0 * pow ( theta, 7 ) / 4725.0;
  
    beta =  2.0                    /     3.0 
          + 2.0 * pow ( theta, 2 ) /    15.0 
          - 4.0 * pow ( theta, 4 ) /   105.0 
          + 2.0 * pow ( theta, 6 ) /   567.0 
          - 4.0 * pow ( theta, 8 ) / 22275.0;

    gamma = 4.0                    /      3.0 
          - 2.0 * pow ( theta, 2 ) /     15.0 
          +       pow ( theta, 4 ) /    210.0 
          -       pow ( theta, 6 ) /  11340.0;
  }
  else
  {
    alpha = ( pow ( theta, 2 ) + theta * sint * cost - 2.0 * sint * sint ) 
      / pow ( theta, 3 );

    beta = ( 2.0 * theta + 2.0 * theta * cost * cost 
      - 4.0 * sint * cost ) / pow ( theta, 3 );

    gamma = 4.0 * ( sint - theta * cost ) / pow ( theta, 3 );
  }

  c2n = + 0.5 * ftab[0] * cos ( t * x[0] );
  for ( i = 2; i < n - 1; i = i + 2 )
  {
    c2n = c2n + ftab[i] * cos ( t * x[i] );
  }
  c2n = c2n + 0.5 * ftab[n-1] * cos ( t * x[n-1] );

  c2nm1 = 0.0;
  for ( i = 1; i <= n - 2; i = i + 2 )
  {
    c2nm1 = c2nm1 + ftab[i] * cos ( t * x[i] );
  }

  value = h * ( 
      alpha * ( ftab[n-1] * sin ( t * x[n-1] )  
              - ftab[0]   * sin ( t * x[0] ) ) 
    + beta * c2n 
    + gamma * c2nm1 );

  free ( x );

  return value;
}
/******************************************************************************/

double filon_fun_sin ( int n, double *f ( int n, double x[] ), double a, 
  double b, double t )

/******************************************************************************/
/*
  Purpose:

    FILON_FUN_SIN uses Filon's method on integrals with a sine factor.

  Discussion:

    The integral to be approximated has the form

      Integral ( A <= X <= B ) F(X) * SIN(T*X) dX

    where T is user specified.

    The function is interpolated over each subinterval by
    a parabolic arc.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 May 2014

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

    Stephen Chase, Lloyd Fosdick,
    An Algorithm for Filon Quadrature,
    Communications of the Association for Computing Machinery,
    Volume 12, Number 8, August 1969, pages 453-457.

    Stephen Chase, Lloyd Fosdick,
    Algorithm 353:
    Filon Quadrature,
    Communications of the Association for Computing Machinery,
    Volume 12, Number 8, August 1969, pages 457-458.

    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.

  Parameters:

    Input, int N, the number of data points, 
    including the endpoints.  N must be odd, and greater than 1.

    Input, external F, the subroutine which evaluates the integrand,
    of the form subroutine F ( N, X, FX ).

    Input, double A, B, the limits of integration.

    Input, double T, multiplier of the X argument of the sine.

    Output, double FILON_FUN_SIN, the approximate value of the integral.
*/
{
  double alpha;
  double beta;
  double cost;
  double *ftab;
  double gamma;
  double h;
  int i;
  double s2n;
  double s2nm1;
  double sint;
  double theta;
  double value;
  double *x;

  if ( a == b )
  {
    value = 0.0;
    return value;
  }
 
  if ( n <= 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "FILON_FUN_SIN - Fatal error!\n" );
    fprintf ( stderr, "  N < 2\n" );
    fprintf ( stderr, "  N = %d\n", n );
    exit ( 1 );
  }
 
  if ( ( n % 2 ) != 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "FILON_FUN_SIN - Fatal error!\n" );
    fprintf ( stderr, "  N must be odd.\n" );
    fprintf ( stderr, "  N = %d\n", n );
    exit ( 1 );
  }
/*
  Set the X values.
*/
  x = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( n - i - 1 ) * a   
           + ( double ) (     i     ) * b ) 
           / ( double ) ( n     - 1 );
  }

  h = ( b - a ) / ( double ) ( n - 1 );
  theta = t * h;
  sint = sin ( theta );
  cost = cos ( theta );

  if ( 6.0 * fabs ( theta ) <= 1.0 )
  {
    alpha = 2.0 * pow ( theta, 3 ) /   45.0 
          - 2.0 * pow ( theta, 5 ) /  315.0 
          + 2.0 * pow ( theta, 7 ) / 4725.0;
  
    beta =  2.0                    /     3.0 
          + 2.0 * pow ( theta, 2 ) /    15.0 
          - 4.0 * pow ( theta, 4 ) /   105.0 
          + 2.0 * pow ( theta, 6 ) /   567.0 
          - 4.0 * pow ( theta, 8 ) / 22275.0;

    gamma = 4.0                    /      3.0 
          - 2.0 * pow ( theta, 2 ) /     15.0 
          +       pow ( theta, 4 ) /    210.0 
          -       pow ( theta, 6 ) /  11340.0;
  }
  else
  {
    alpha = ( pow ( theta, 2 ) + theta * sint * cost 
      - 2.0 * sint * sint ) / pow ( theta, 3 );

    beta = ( 2.0 * theta + 2.0 * theta * cost * cost
      - 4.0 * sint * cost ) / pow ( theta, 3 );

    gamma = 4.0 * ( sint - theta * cost ) / pow ( theta, 3 );
  }
/*
  Tabulate the function.
*/
  ftab = f ( n, x );

  s2n = + 0.5 * ftab[0] * sin ( t * x[0] );
  for ( i = 2; i < n - 1; i = i + 2 )
  {
    s2n = s2n + ftab[i] * sin ( t * x[i] );
  }
  s2n = s2n + 0.5 * ftab[n-1] * sin ( t * x[n-1] );

  s2nm1 = 0.0;
  for ( i = 1; i <= n - 2; i = i + 2 )
  {
    s2nm1 = s2nm1 + ftab[i] * sin ( t * x[i] );
  }

  value = h * ( 
      alpha * ( ftab[0]   * cos ( t * x[0] ) 
              - ftab[n-1] * cos ( t * x[n-1] ) )
    + beta * s2n 
    + gamma * s2nm1 );
 
  free ( ftab );
  free ( x );

  return value;
}
/******************************************************************************/

double filon_tab_sin ( int n, double ftab[], double a, double b, double t )

/******************************************************************************/
/*
  Purpose:

    FILON_TAB_SIN uses Filon's method on integrals with a sine factor.

  Discussion:

    The integral to be approximated has the form

      Integral ( A <= X <= B ) F(X) * SIN(T*X) dX

    where T is user specified.

    The function is interpolated over each subinterval by
    a parabolic arc.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 May 2014

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

    Stephen Chase, Lloyd Fosdick,
    An Algorithm for Filon Quadrature,
    Communications of the Association for Computing Machinery,
    Volume 12, Number 8, August 1969, pages 453-457.

    Stephen Chase, Lloyd Fosdick,
    Algorithm 353:
    Filon Quadrature,
    Communications of the Association for Computing Machinery,
    Volume 12, Number 8, August 1969, pages 457-458.

    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.

  Parameters:

    Input, int N, the number of data points, 
    including the endpoints.  N must be odd, and greater than 1.

    Input, double FTAB[N], contains the value of the function
    at A, A+H, A+2*H, ... , B-H, B, where H = (B-A)/(N-1).

    Input, double A, B, the limits of integration.

    Input, double T, multiplier of the X argument of the sine.

    Output, double FILON_TAB_SIN, the approximate value of the integral.
*/
{
  double alpha;
  double beta;
  double cost;
  double gamma;
  double h;
  int i;
  double s2n;
  double s2nm1;
  double sint;
  double theta;
  double value;
  double *x;

  if ( a == b )
  {
    value = 0.0;
    return value;
  }
 
  if ( n <= 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "FILON_TAB_SIN - Fatal error!\n" );
    fprintf ( stderr, "  N < 2\n" );
    fprintf ( stderr, "  N = %d\n", n );
    exit ( 1 );
  }
 
  if ( ( n % 2 ) != 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "FILON_TAB_SIN - Fatal error!\n" );
    fprintf ( stderr, "  N must be odd.\n" );
    fprintf ( stderr, "  N = %d\n", n );
    exit ( 1 );
  }
/*
  Set the X values.
*/
  x = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( n - i - 1 ) * a   
           + ( double ) (     i     ) * b ) 
           / ( double ) ( n     - 1 );
  }

  h = ( b - a ) / ( double ) ( n - 1 );
  theta = t * h;
  sint = sin ( theta );
  cost = cos ( theta );

  if ( 6.0 * fabs ( theta ) <= 1.0 )
  {
    alpha = 2.0 * pow ( theta, 3 ) /   45.0 
          - 2.0 * pow ( theta, 5 ) /  315.0 
          + 2.0 * pow ( theta, 7 ) / 4725.0;
  
    beta =  2.0                    /     3.0 
          + 2.0 * pow ( theta, 2 ) /    15.0 
          - 4.0 * pow ( theta, 4 ) /   105.0 
          + 2.0 * pow ( theta, 6 ) /   567.0 
          - 4.0 * pow ( theta, 8 ) / 22275.0;

    gamma = 4.0                    /      3.0 
          - 2.0 * pow ( theta, 2 ) /     15.0 
          +       pow ( theta, 4 ) /    210.0 
          -       pow ( theta, 6 ) /  11340.0;
  }
  else
  {
    alpha = ( pow ( theta, 2 ) + theta * sint * cost 
      - 2.0 * sint * sint ) / pow ( theta, 3 );

    beta = ( 2.0 * theta + 2.0 * theta * cost * cost 
      - 4.0 * sint * cost ) / pow ( theta, 3 );

    gamma = 4.0 * ( sint - theta * cost ) / pow ( theta, 3 );
  }
  
  s2n = + 0.5 * ftab[0] * sin ( t * x[0] );
  for ( i = 2; i < n - 1; i = i + 2 )
  {
    s2n = s2n + ftab[i] * sin ( t * x[i] );
  }
  s2n = s2n + 0.5 * ftab[n-1] * sin ( t * x[n-1] );

  s2nm1 = 0.0;
  for ( i = 1; i <= n - 2; i = i + 2 )
  {
    s2nm1 = s2nm1 + ftab[i] * sin ( t * x[i] );
  }

  value = h * ( 
      alpha * ( ftab[0]   * cos ( t * x[0] ) 
              - ftab[n-1] * cos ( t * x[n-1] ) ) 
    + beta * s2n 
    + gamma * s2nm1 );
 
  free ( x );

  return value;
}
/******************************************************************************/

int i4_power ( int i, int j )

/******************************************************************************/
/*
  Purpose:

    I4_POWER returns the value of I^J.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 October 2007

  Author:

    John Burkardt

  Parameters:

    Input, int I, J, the base and the power.  J should be nonnegative.

    Output, int I4_POWER, the value of I^J.
*/
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "I4_POWER - Fatal error!\n" );
      fprintf ( stderr, "  I^J requested, with I = 0 and J negative.\n" );
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "I4_POWER - Fatal error!\n" );
      fprintf ( stderr, "  I^J requested, with I = 0 and J = 0.\n" );
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
}
/******************************************************************************/

void timestamp ( )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}


double rand0_1()
{
  //Return a random number between the limits {0,1}
  double rnd;
  rnd = (double)rand() / double(RAND_MAX);
  return rnd;
}

double rand_unif_rms( double rms )
{
  //Return a random number between the limits {-sqrt(3)rms, +sqrt(3)rms}
  // with uniform distribution with having root-mean square = rms 
  double rnd;
  rnd = (double)rand() / double(RAND_MAX);
  rnd = 2.*sqrt(3.)*rms*(rnd-0.5);
  return rnd;
}

double rand_gaus_rms( double rms )
{
  // Gaussian distribution with rms
  double rnd1, rnd2;
  rnd1 = (double)rand() / double(RAND_MAX);
  rnd2 = (double)rand() / double(RAND_MAX);
  return rms*sqrt(2.*log(1/rnd1))*cos(6.283185307179586*rnd2);
}

double rand_gaus_rms_trunc( double rms, double rms_lim )
{
  // Gaussian distribution with rms and truncated at |rms_lim|
  double rnd1, rnd2, val=(rms_lim+1.);
  while( fabs(val)>rms_lim ){
    rnd1 = (double)rand() / double(RAND_MAX);
    rnd2 = (double)rand() / double(RAND_MAX);
    val = rms*sqrt(2.*log(1/rnd1))*cos(6.283185307179586*rnd2);
  }
  return val;
}


/******************************************************************************/
double *FourierTransCos( int Nk, double k[], int Nx, double x[], double fx[])
{
  double *FTf;
  int j;
  FTf = (double*) malloc(Nk*sizeof(double));
  for ( j=0; j<Nk ; j++ ) {
    FTf[j]=0.39894228*filon_tab_cos(Nx,fx,x[0],x[Nx-1],k[j]);
  }
  return FTf ;
}
double *FourierTransSin( int Nk, double k[], int Nx, double x[], double fx[])
{
  double *FTf;
  int j;
  FTf = (double*) malloc(Nk*sizeof(double));
  for ( j=0; j<Nk ; j++ ) {
    FTf[j]=0.39894228*filon_tab_sin(Nx,fx,x[0],x[Nx-1],k[j]);
  }
  return FTf ;
}
double *InvFourierTransRe( int Nk, double k[], double FTfRe[], double FTfIm[],
  int Nx, double x[])
{
  double *IFTfRe;
  int i;
  IFTfRe = (double*) malloc(Nx*sizeof(double));
  for ( i=0; i<Nx ; i++ ) {
    IFTfRe[i]= filon_tab_cos(Nk,FTfRe,k[0],k[Nk-1],x[i]);
    IFTfRe[i]+=filon_tab_sin(Nk,FTfIm,k[0],k[Nk-1],x[i]);
    IFTfRe[i]*=0.39894228;
  }
  return IFTfRe ;
}
double *InvFourierTransIm( int Nk, double k[], double FTfRe[], double FTfIm[],
  int Nx, double x[])
{
  double *IFTfIm;
  int i;
  IFTfIm = (double*) malloc(Nx*sizeof(double));
  for ( i=0; i<Nx ; i++ ) {
    IFTfIm[i]= filon_tab_sin(Nk,FTfRe,k[0],k[Nk-1],x[i]);
    IFTfIm[i]-=filon_tab_cos(Nk,FTfIm,k[0],k[Nk-1],x[i]);
    IFTfIm[i]*=0.39894228;
  }
  return IFTfIm ;
}

double *p_factor_radial( int Nk, double k[], int n, int m_max, 
    double fk[], double r, double r0 )
  {
    double *FTb_n_0;
    double p_m, q_n_m, sum;
    int j, m, mm;
    if( n==0 ) {printf("n=0 not allowed in p_factor_radial\n"); exit(-1); }
    FTb_n_0 = (double*) malloc(Nk*sizeof(double));
    for ( j=0; j<Nk ; j++ ) {
      sum = 0;
      for( m=0; m<m_max; m++){
        p_m = 1 ;
        if( m>0 ){
          for( mm=1; mm<=m ; mm++) {
            q_n_m = (double) (n+2*mm)/(4*mm*(n+mm))/(n+2*mm-2);
            q_n_m *= (r0*r0*k[j]*k[j]) ;
            p_m *= q_n_m ;
          }
        }
        //sum += p_m*pow( (r/r0), n+2*m-1) ;
        sum += p_m*pow( (r/r0), n+2*m-2) ; //To obtain correct scaling
        // Correct scaling refers to obtaining the same b_n,0 result,
        // independent of r (position where fields are measured).
      }
      FTb_n_0[j] = fk[j]/sum;
    }
    //printf(" r/r0 = %lf \n", r/r0); exit(0);
    return FTb_n_0 ;
  }

double *p_factor_tang( int Nk, double k[], int n, int m_max, 
    double fk[], double r, double r0 )
  {
    double *b_n_0;
    double p_m, q_n_m, sum;
    int j, m, mm;
    b_n_0 = (double*) malloc(Nk*sizeof(double));
    for ( j=0; j<Nk ; j++ ) {
      sum = 0;
      if( n==0 ) {printf("n=0 not allowed in p_factor_radial\n"); exit(-1); }
      for( m=0; m<m_max; m++){
        p_m = 1 ;
        if( m>0 ){
          for( mm=1; mm<(m+1) ; mm++) {
            q_n_m = (double) (n+2*mm)/(4*mm*(n+mm))/(n+2*mm-2);
            q_n_m *= (r0*r0*k[j]*k[j]) ;
            p_m *= q_n_m ;
          }
        }
        sum += p_m * (double) ( (n)/(n+2*m) ) ;
      }
      b_n_0[j] = fk[j]/sum;
    }
    return b_n_0 ;
  }


double integd_trapez( int n, double x[], double ftab[] )
{
  //Simple trapezoidal integration of ftab(x) for two arrays as input
  int i;
  double value = 0.0;
  if ( n < 2 ) { printf(" trapzd ERROR. n<2 \n"); exit(-1); }
  if ( x[0] == x[1] ) return value;
  if ( n == 2 ) {
    value = 0.5*( ftab[1]+ftab[0] )*( x[1] - x[0] );
  } else {
    for ( i=1; i<n ; i++) value += 0.5*( ftab[i]+ftab[i-1] )*( x[i] - x[i-1] );
  }
  return value;
}

double *derivative_primative_d( int n, double x[], double fx[])
{
  int i;
  double *dfx;
  dfx = (double*) malloc(n*sizeof(double));
  dfx[0] = (fx[1] - fx[0])/(x[1] - x[0]);
  dfx[n-1] = (fx[n-1] - fx[n-2])/(x[n-1] - x[n-2]);
  for( i=1; i<(n-1); i++)
  {
    dfx[i] = (fx[i+1] - fx[i-1])/(x[i+1] - x[i-1]);
  }
  return dfx ;
}

double *derivative_d( int n, double x[], double fx[])
{
  //Derivation that applied Richardson's extrapolation for error in O(h^4)
  // Boundaries use more conventional methods where error are more significant.
  int i;
  double *dfx;
  dfx = (double*) malloc(n*sizeof(double));
  if ( n==1 ) { printf(" derivative_d: Error, n<= 1.\n"); exit(-1); }
  dfx[0] = (fx[1] - fx[0])/(x[1] - x[0]);
  dfx[n-1] = (fx[n-1] - fx[n-2])/(x[n-1] - x[n-2]);
  if ( n<=2 ) return dfx ;
  i=1;
  dfx[i] = (fx[i+1] - fx[i-1])/(x[i+1] - x[i-1]);
  i=(n-2);
  dfx[i] = (fx[i+1] - fx[i-1])/(x[i+1] - x[i-1]);
  for( i=2; i<(n-2); i++)
  {
    dfx[i] = -fx[i+2] + 8*fx[i+1] - 8*fx[i-1] + fx[i-2] ;
    dfx[i] /= 3*(x[i+2] - x[i-2]);
  }
  return dfx ;
}