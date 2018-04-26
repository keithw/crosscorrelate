/*
   program to calculate cross-correlation between two mahimahi trace files

   Usage: ./xcor BIN_DURATION (in milliseconds) file1 file2

   (Can set file1=file2 to get autocorrelation.)

   xcor will first convert mahimahi's traces (which record the arrival
   time of each packet, in milliseconds) to "throughput" in packets
   per bin. E.g. if bin = 100, xcor will calculate throughput over
   time in packets per tenth-of-a-second interval.

   Next, xcor calculates the cross-correlation (at lags of -1 minute to +1 minute)
   between the two binned throughput traces. This is printed on
   standard output.

   compile with: g++ -std=c++17 -g -O2 -Wpedantic -Wall -Wextra -Weffc++ -Werror -o xcor xcor.cc

   (requires a C++17 compiler, such as recent GCC 7 or recent clang)

   keithw@cs.stanford.edu 2018-4-25
*/

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <numeric>
#include <cmath>
#include <fstream>

using namespace std;

/* parse string representing an unsigned integer and verify that it round-trips */
int careful_atoi( const string & str );

/* read a sequence of integers from an input stream */
vector<int> read_integer_sequence( istream & in );

/* bin a sequence of integers into buckets, each bin_duration in width */
vector<int> aggregate( const vector<int> & events, const int bin_duration );

/* calculate mean and variance of a sequence of integers */
pair<double, double> statistics( const vector<int> & input );

/* cross-correlate two sequence of integers up to some max lag */
vector<pair<int, double>> crosscorrelate( const vector<int> & input1,
					  const vector<int> & input2,
					  const int max_lag );

/* open a file and throw exception if fails */
class careful_ifstream : public ifstream
{
public:
  careful_ifstream( const string & filename ) : ifstream( filename )
  { if ( not good() ) { throw runtime_error( "can't open " + filename ); } }
};

int main( int argc, char *argv[] )
{
  if ( argc <= 0 ) {
    abort();
  } else if ( argc != 4 ) {
    cerr << "Usage: " << argv[ 0 ] << " BIN_DURATION (in milliseconds) trace1 trace2\n";
    return EXIT_FAILURE;
  }

  try {
    const int bin_duration = careful_atoi( argv[ 1 ] );

    /* open trace files */
    careful_ifstream trace1 { argv[ 2 ] }, trace2 { argv[ 3 ] };

    /* read trace files into memory */
    const auto arrival_times1_ms = read_integer_sequence( trace1 );
    const auto arrival_times2_ms = read_integer_sequence( trace2 );

    /* convert packet arrival times into throughput bins */
    const auto arrival_counts1 = aggregate( arrival_times1_ms, bin_duration );
    const auto arrival_counts2 = aggregate( arrival_times2_ms, bin_duration );

    /* cross-correlate the throughput bins */
    const auto crosscorrelation = crosscorrelate( arrival_counts1,
						  arrival_counts2,
						  60000 / bin_duration ); /* +/- one second of lag */

    /* print out the cross-correlation as a function of lag */
    for ( const auto & x : crosscorrelation ) {
      cout << x.first * static_cast<signed>( bin_duration ) << ": " << x.second << "\n";
    }
  } catch ( const exception & e ) {
    cerr << argv[ 0 ] << ": " << e.what() << "\n";
    return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;
}

int careful_atoi( const string & str )
{
  const int ret = stoul( str );

  if ( str != to_string( ret ) ) {
    throw runtime_error( "invalid int: " + str );
  }

  return ret;
}

vector<int> read_integer_sequence( istream & in )
{
  vector<int> ret;
  
  /* read the trace from stdin */
  string line;
  while ( true ) {
    getline( in, line );
    if ( not cin.good() or line.empty() ) { break; }

    ret.emplace_back( careful_atoi( line ) );
  }

  return ret;
}

vector<int> aggregate( const vector<int> & events, const int bin_duration )
{
  if ( events.empty() ) { throw runtime_error( "can't bin empty list of events" ); }

  vector<int> ret( 1 + events.back() / bin_duration );

  for ( const int event_time : events ) {
    ret.at( event_time / bin_duration )++;
  }
  
  return ret;
}

pair<double, double> statistics( const vector<int> & input )
{
  if ( input.empty() ) { throw runtime_error( "can't calculate statistics on empty vector" ); }

  const double sum = accumulate( input.begin(), input.end(), 0 );
  const double mean = sum / input.size();

  double total_difference = 0;
  double total_variance = 0;
  for ( const auto value : input ) {
    const double diff = value - mean;
    total_difference += diff;
    total_variance += diff * diff;
  }

  const double variance = (total_variance - total_difference * total_difference / input.size()) / (input.size() - 1);

  return { mean, variance };
}

vector<pair<int, double>> crosscorrelate( const vector<int> & input1,
					  const vector<int> & input2,
					  const int max_lag )
{
  vector<pair<int, double>> ret( 2 * max_lag + 1 );

  const auto [ mean1, variance1 ] = statistics( input1 );
  const auto [ mean2, variance2 ] = statistics( input2 );

  for ( int lag = -max_lag; lag <= max_lag; lag++ ) {
    auto & output = ret.at( lag + max_lag );
    output.first = lag;
    
    double sum = 0;
    int count = 0;
    for ( int index1 = 0; index1 < signed( input1.size() ); index1++ ) {
      const int index2 = index1 + lag;
      if ( index2 >= 0 and index2 < signed( input2.size() ) ) {
	sum += (input1.at( index1 ) - mean1) * (input2.at( index2 ) - mean2);
	count++;
      }
    }

    output.second = (sum / count) / ( sqrt( variance1 ) * sqrt( variance2 ) );
  }
  
  return ret;
}
