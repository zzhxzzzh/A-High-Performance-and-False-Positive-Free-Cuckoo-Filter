#ifndef HCF_HCF
#define HCF_HCF



#include <thread>  
#include <atomic>
#include <assert.h>
#include <algorithm>
#include <thread>
#include <atomic>
#include <future> 

#include "hashutil.h"
#include "singletable.h"
#include <bitset>
#include <iostream>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <functional>
#include <cstring>
#include <smmintrin.h>
#include "crc.h"
#include <cmath> 
#include <iostream>
#include <immintrin.h>  

using namespace std;

namespace HCF {

template <typename T>  
uint32_t CRC32Hash1(T key, uint32_t seed) 
{

    const char* data = reinterpret_cast<const char*>(&key);
    size_t length = sizeof(T);
    uint32_t crc = seed;

    while (length >= 4) {
        crc = _mm_crc32_u32(crc, *reinterpret_cast<const uint32_t*>(data));
        data += 4;
        length -= 4;
    }

    while (length--) {
        crc = _mm_crc32_u8(crc, *data++);
    }

    return crc;

}


template <typename T>  
uint32_t CRC32Hash2(T key, uint32_t seed) 
{
    char* k= (char*) malloc(sizeof(T));
    k = (char*) memcpy(k,&key,sizeof(T));

    uint32_t r=crc32(k,sizeof(T),seed);
    free(k);
    return (uint32_t) r;

}

enum Status {
  Ok = 0,              
  NotFound = 1,        
  NotEnoughSpace = 2,  
  NotSupported = 3,    
};


const size_t kMaxCuckooCount = 500;


template <typename ItemType, size_t bits_per_item,template <size_t> class TableType = SingleTable,typename HashFamily = TwoIndependentMultiplyShift>
class HCF1 {
  TableType<bits_per_item> *table_;
  size_t num_items_;
  int fp_size;
  int p;
  int q; 

  int distance;

  uint32_t dis_bit=ceil(log2(p));

 
  typedef struct {
    size_t index;   
    uint32_t distance;   
    uint32_t fp;    
    bool used;     
  } VictimCache;

  VictimCache victim_;  

  HashFamily hasher_;  

 public:


 
  inline size_t IndexHash(uint32_t hv) const {
    hv=(hv>>fp_size); 
    return hv %(table_->NumBuckets());  
  }


  inline uint32_t FpHash(uint32_t hv) const {
    uint32_t mask = (1U << fp_size) - 1; 
    return  ((hv & mask) == 0 ? 1 : (hv & mask)); 
  }


 
  inline void GenerateIndexFpHash(const ItemType& item, size_t* index,uint32_t* fp) const {
    uint32_t hash_value =CRC32Hash1<ItemType>(item,444);   
     *index = IndexHash(hash_value);  
     *fp = FpHash(hash_value);           
  }


  inline size_t AltIndex(const size_t index, const uint32_t fp) const {
    int hash_value=CRC32Hash2<uint32_t>(fp,444);
    return ((uint32_t)(index ^ (hash_value)))%(table_->NumBuckets());
  }


  inline size_t final_num_buckets()const{
    return table_->NumBuckets();
  }


  Status AddImpl(const size_t i, const uint32_t fp);

  Status AddImpl2(const size_t i, const uint32_t fp);



  inline size_t upperpower2(size_t num) {
    size_t power = 1;
    while (power < num) {
        power <<= 1;
    }
    return power;
  }



  double LoadFactor() const { return 1.0 * Size() / table_->SizeInTags(); }


 
  double BitsPerItem() const { return 8.0 * table_->SizeInBytes() / Size(); }



  explicit HCF1(const size_t max_num_keys) : num_items_(0), victim_(), hasher_(),fp_size(8),p(16),q(12){
    size_t assoc = 4;  
    size_t num_buckets = upperpower2(std::max<uint64_t>(1, max_num_keys / assoc));  
    double frac = (double)max_num_keys / num_buckets / assoc;

    if (frac > 0.96) {
      num_buckets <<= 1;  
    }
    victim_.used = false;
    table_ = new TableType<bits_per_item>(num_buckets);  

  }


  
  Status Add(const ItemType &item);

  
  ~HCF1 () { delete table_; }

 

  Status Contain1(const ItemType &item) const;

  Status Contain2(const ItemType &item,const std::vector<std::vector<int32_t>>& all_distances) const;

  
  Status ContainWithIndexFp(size_t i1, uint32_t fp) const;

  
  Status Delete(const ItemType &item);

 
  Status Deletewithindexfp(size_t i1, uint32_t fp);


  std::string Info() const;


  size_t Size() const { return num_items_; }


  size_t SizeInBytes() const { return table_->SizeInBytes(); }

  vector<std::vector<int32_t>> get_ReadAllDistances(const ItemType &key){
    bool found = false;
    size_t i1, i2;
    uint32_t fp;

    GenerateIndexFpHash(key, &i1, &fp);  
    i2= AltIndex(i1, fp);
  
    return table_-> ReadAllDistances(i1,i2,q,p);
  }


  
  void PrintAllTags() const {
    table_->PrintAllTags1();
    
  }

  void get_ReadFingerprint(const size_t i,const size_t j) const{
    table_->testReadFingerprint(i,j);
  }



};



template <typename ItemType, size_t bits_per_item,template <size_t> class TableType, typename HashFamily>
Status HCF1 <ItemType, bits_per_item, TableType, HashFamily>::Add(const ItemType &item) {

  size_t i;
  uint32_t fp;
  if (victim_.used) {
    return NotEnoughSpace;  
  }
 
  GenerateIndexFpHash(item, &i, &fp);  
  return AddImpl(i, fp); 
}


template <typename ItemType, size_t bits_per_item, template <size_t> class TableType, typename HashFamily>
Status HCF1<ItemType, bits_per_item, TableType, HashFamily>::AddImpl(const size_t i,const uint32_t fp) {
    size_t curindex = i;  
    size_t curfp=fp;
    uint32_t oldfp; 
    uint32_t oldindex; 
    uint32_t olddistance; 


    
    for (uint32_t count = 0; count < kMaxCuckooCount; count++) {
        
      bool kickout = (count) > 0; 
      oldfp = 0;
      
      
      for (size_t attempt = 0; attempt < q; ++attempt) {
          size_t new_bucket_i = (curindex + attempt) % table_->num_buckets_; 
          for (size_t j = 0; j < table_->kTagsPerBucket; j++) {
              if (table_->ReadFingerprint(new_bucket_i, j) == 0 && table_->ReadDistance(new_bucket_i, j)==0) {
                  table_->WriteFingerprint(new_bucket_i, j, curfp); 
                  table_->WriteDistance(new_bucket_i, j,attempt); 
                  num_items_++; 
                  return Ok;
              }
          }
      }


      size_t A_curindex = AltIndex(curindex, curfp); 
      for (size_t attempt = q; attempt < p; ++attempt) {
          size_t new_bucket_i = (A_curindex + attempt - q) % table_->num_buckets_;
          
          for (size_t j = 0; j < table_->kTagsPerBucket; j++) {
              if (table_->ReadFingerprint(new_bucket_i, j) == 0 && table_->ReadDistance(new_bucket_i, j)==0) {
                  table_->WriteFingerprint(new_bucket_i, j, curfp); 
                  table_->WriteDistance(new_bucket_i, j,attempt); 
                  num_items_++; 
                  return Ok;
              }
          }
      }
      


      
      if (kickout) {
          if (rand() % 2 == 0) {
            curindex = AltIndex(curindex,curfp);
            size_t kickout_i=(curindex+(rand() % (p - q)))%table_->num_buckets_;
            size_t r = rand() % table_->kTagsPerBucket;
            oldfp = table_->ReadFingerprint(kickout_i, r);
            olddistance=table_->ReadDistance(kickout_i, r);
            oldindex = kickout_i;

           
            table_->WriteFingerprint(kickout_i, r,curfp); 
            table_->WriteDistance(kickout_i,r,kickout_i-curindex+q);

            curfp = oldfp; 
            
            if(olddistance<q){
                curindex = oldindex-olddistance;
            }else{
              curindex =AltIndex((oldindex-(olddistance-q)),curfp);
            }       
          }else {
           
            size_t kickout_i=(curindex+(rand()%q))%table_->num_buckets_;
            size_t r = rand() % table_->kTagsPerBucket;

            oldfp = table_->ReadFingerprint(kickout_i, r); 
            olddistance=table_->ReadDistance(kickout_i, r);
            oldindex = kickout_i;

           
            table_->WriteFingerprint(kickout_i, r,curfp);
            table_->WriteDistance(kickout_i,r,kickout_i-curindex);


            curfp = oldfp;

            if(olddistance>=q){
                curindex =AltIndex((oldindex-(olddistance-q)),curfp);
            }else{
                curindex = oldindex-olddistance;
            }
          }
      }

    }

    victim_.index = curindex;
    victim_.fp =curfp;
    victim_.distance=olddistance;
    victim_.used = true;
    return Ok;
}




template <typename ItemType, size_t bits_per_item,template <size_t> class TableType, typename HashFamily>
Status HCF1<ItemType, bits_per_item, TableType, HashFamily>::Contain1(const ItemType &key) const {
  bool found = false;
  size_t i1, i2;
  uint32_t fp;

  GenerateIndexFpHash(key, &i1, &fp);  
  

  found = victim_.used && (fp == victim_.fp) &&(i1 == victim_.index || i2 == victim_.index);


  if (found || table_-> Find_fp_In_Main_Buckets(i1,fp,q,p)) {
    return Ok;  
  }else{
    i2 = AltIndex(i1, fp);  

    if(table_-> Find_fp_In_Auxiliary_Buckets(i2,fp,q,p)){
      return Ok;
    }else{
       return NotFound;  
    }
  }
  
}


template <typename ItemType, size_t bits_per_item, template <size_t> class TableType, typename HashFamily>
Status HCF1<ItemType, bits_per_item, TableType, HashFamily>::ContainWithIndexFp(size_t i1, uint32_t fp) const {
    bool found = false;
    size_t i2;  


    found = victim_.used && (fp == victim_.fp) && (i1 == victim_.index || i2 == victim_.index);

    if (found || table_->Find_fp_In_Main_Buckets(i1, fp, q, p)) {
        return Ok;  
    } else {
        i2 = AltIndex(i1, fp); 
        if(table_-> Find_fp_In_Auxiliary_Buckets(i2,fp,q,p)){
          return Ok;
        }else{
          return NotFound;  
        }
    }
}





template <typename ItemType, size_t bits_per_item,template <size_t> class TableType, typename HashFamily>
Status HCF1<ItemType, bits_per_item, TableType, HashFamily>::Delete(const ItemType &key) {
  size_t i1, i2;
  uint32_t fp;

  GenerateIndexFpHash(key, &i1, &fp);  
  

  if (table_->DeleteTagFrom_MainBuckets(i1, fp,q,p)) {
    num_items_--;  
    goto TryEliminateVictim;
  } else{
      i2 = AltIndex(i1, fp);  
      if (table_->DeleteTagFrom_AuxiliaryBuckets(i2, fp,q,p)) {
          num_items_--;
          goto TryEliminateVictim;
      } else if (victim_.used && fp == victim_.fp &&(i1 == victim_.index || i2 == victim_.index)) {
            
            victim_.used = false;
            return Ok;
      } else {
          return NotFound;  
      }
  }

  TryEliminateVictim:
    if (victim_.used) {
        victim_.used = false;
        size_t i = victim_.index;
        uint32_t fp = victim_.fp;
        AddImpl(i, fp);  
    }
    return Ok;
 
}



template <typename ItemType, size_t bits_per_item,template <size_t> class TableType, typename HashFamily>
Status HCF1<ItemType, bits_per_item, TableType, HashFamily>::Deletewithindexfp(size_t i1, uint32_t fp) {
  size_t  i2;

  if (table_->DeleteTagFrom_MainBuckets(i1, fp,q,p)) {
    num_items_--; 
    goto TryEliminateVictim;
  } else{
      i2 = AltIndex(i1, fp);  
      if (table_->DeleteTagFrom_AuxiliaryBuckets(i2, fp,q,p)) {
          num_items_--;
          goto TryEliminateVictim;
      } else if (victim_.used && fp == victim_.fp &&(i1 == victim_.index || i2 == victim_.index)) {
            
            victim_.used = false;
            return Ok;
      } else {
          return NotFound; 
      }
  }

  TryEliminateVictim:
    if (victim_.used) {
        victim_.used = false;
        size_t i = victim_.index;
        uint32_t fp = victim_.fp;
        AddImpl(i, fp);  
    }
    return Ok;
 
}


template <typename ItemType, size_t bits_per_item, template <size_t> class TableType, typename HashFamily>
std::string HCF1<ItemType, bits_per_item, TableType, HashFamily>::Info() const {
  std::stringstream ss;
  ss << "HCF Status:\n"
     << "\t\t" << table_->Info() << "\n"
     << "\t\tCurrent number of inserted elements: " << Size() << "\n"
     << "\t\tLoad factor: " << LoadFactor() << "\n"
     << "\t\tTotal table size in bytes: " << (table_->SizeInBytes()) << " B\n";
  if (Size() > 0) {
    ss << "\t\tBits per item: " << BitsPerItem() << "\n";
  } else {
    ss << "\t\tBits per item: N/A\n";
  }
  return ss.str();
}


}  // namespace HCF

#endif  // HCF
