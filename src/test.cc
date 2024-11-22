#include "hcf.h"  

#include <assert.h>        
#include <math.h>          

#include <iostream>       
#include <vector>          
#include "random.h"

using namespace std;

using HCF::HCF1;  

int main(int argc, char **argv) {
  size_t total_items = 1000000;  
  HCF1<size_t, 8> filter(total_items);

  // Generate random input data
  const vector<uint32_t> input = GenerateRandom32(2 * total_items);

  // Insert 
  for (size_t i = 0; i < total_items; ++i) {
      if (filter.Add(input[i]) != HCF::Ok) {
          cout << "Insertion failed: (" << input[i] << ", " << i << ")" << endl;
          break;
      }
  }


  // Verify all inserted items can be found
  for (size_t i = 0; i < total_items; ++i) {
      assert(filter.Contain1(input[i]) == HCF::Ok);
  }

  cout << "Insertion and verification completed!\n";

  // Delete some items
  for (size_t i = 0; i < total_items / 10; ++i) {
      if (filter.Delete(input[i]) == HCF::Ok) {
          //cout << "Deletion successful: (" << input[i] << ", " << i << ")" << endl;
      } else {
          //cout << "Deletion failed: (" << input[i] << ", " << i << ")" << endl;
      }
  }

  // Verify that the deleted items are no longer present
  for (size_t i = 0; i < total_items / 10; ++i) {
      if (filter.Contain1(input[i]) == HCF::Ok) {
          cout << "Deleted item found: (" << input[i] << ", " << i << ")" << endl;
      }
  }


  // Test false positive rate
  size_t false_queries = 0;
  size_t total_queries = 0;
  for (size_t i = total_items; i < 2 * total_items; ++i) {
      if (filter.Contain1(input[i]) == HCF::Ok) {
          false_queries++;
      }
      total_queries++;
  }

  cout << "Total queries: " << total_queries << endl;
  cout << "False positive rate: " << (100.0 * false_queries / total_queries) << "%\n";

  cout << "Program executed successfully!\n";


  return 0;
}
