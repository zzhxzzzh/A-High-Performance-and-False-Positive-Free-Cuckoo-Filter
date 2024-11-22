#ifndef HCF_HNFP
#define HCF_HNFP

#include <assert.h>

#include <sstream>
#include <cstring> 
#include <utility>  
#include <cstdint> 
#include <iostream>
#include <iomanip>
#include <cstdint>
#include <utility> 
#include <x86intrin.h> 
#include <immintrin.h> 
#include <chrono>
#include <emmintrin.h> 
#include <smmintrin.h> 
#include <thread>
#include <vector>
#include <iostream>
#include <atomic>


using namespace std;

namespace HCF {


template <size_t bits_per_fp>
class SingleTable {

 public:  
 
  static const size_t kTagsPerBucket = 4;

  static const size_t bits_per_distance=4;

  static const size_t kBytesPerBucket =
      (bits_per_fp * kTagsPerBucket + 7) >> 3;

  

  static const size_t kBytesPerBucketDistance =
      (bits_per_distance * kTagsPerBucket + 7) >> 3;

  static const uint32_t kFpMask = (1 << bits_per_fp) - 1;
  static const uint32_t kDistanceMask = (1 << bits_per_distance) - 1;


  static const size_t kPaddingBuckets =
    ((((kBytesPerBucket + 7) / 8) * 8) - 1) / kBytesPerBucket;

  static const size_t kPaddingBucketsDistance =
    ((((kBytesPerBucketDistance + 7) / 8) * 8) - 1) / kBytesPerBucketDistance;  


  struct Bucket {
    char bits_[kBytesPerBucket];
  } __attribute__((__packed__)); 



  struct BucketDistance {
    char bits_[kBytesPerBucketDistance];
  } __attribute__((__packed__)); 



  Bucket *buckets_fp_;


  BucketDistance *buckets_distance_ ;


 
  size_t num_buckets_;



  explicit SingleTable(const size_t num): num_buckets_(num) {
       
        buckets_fp_ = new Bucket[num_buckets_ + kPaddingBuckets];
      
        memset(buckets_fp_, 0, kBytesPerBucket * (num_buckets_ + kPaddingBuckets));

        
        buckets_distance_ = new BucketDistance[num_buckets_ + kPaddingBucketsDistance];
       
        memset(buckets_distance_,0, kBytesPerBucketDistance * (num_buckets_ + kPaddingBucketsDistance));
    }



  ~SingleTable() { 
    delete[] buckets_fp_;
    delete[] buckets_distance_;
  }


  size_t NumBuckets() const {
    return num_buckets_;
  }


  size_t SizeInBytes() const { 
    return kBytesPerBucket * num_buckets_; 
  }


  size_t SizeInTags() const { 
    return kTagsPerBucket * num_buckets_; 
  }





std::string Info() const {
    std::stringstream ss;
    ss << "Bits per fingerprint in the table: " << bits_per_fp << " bits\n";
    ss << "\t\tNumber of slots per bucket: " << kTagsPerBucket << "\n";
    ss << "\t\tNumber of buckets: " << num_buckets_ << "\n";
    ss << "\t\tTotal number of slots: " << SizeInTags() << "\n";
    return ss.str();
}


  inline uint32_t ReadFingerprint(const size_t i, const size_t j) const {
    const char *p = buckets_fp_[i].bits_; 
    uint32_t fp; 

    if (bits_per_fp == 2) {
      fp = *((uint8_t *)p) >> (j * 2);
    } else if (bits_per_fp == 4) {
      p += (j >> 1); 
      fp = *((uint8_t *)p) >> ((j & 1) << 2);
    } else if (bits_per_fp == 8) {
      p += j; 
      fp = *((uint8_t *)p);
    } else if (bits_per_fp == 12) {
      p += j + (j >> 1); 
      fp = *((uint16_t *)p) >> ((j & 1) << 2);
    } else if (bits_per_fp == 16) {
      p += (j << 1); 
      fp = *((uint16_t *)p);
    } else if (bits_per_fp == 32) {
      fp = ((uint32_t *)p)[j];
    }
    return fp & kFpMask; 
  }


  inline uint32_t ReadDistance(const size_t i, const size_t j) const {
 
    const char *p = buckets_distance_[i].bits_;
    uint32_t distance; 

    if (bits_per_distance == 2) {
      distance = *((uint8_t *)p) >> (j * 2);
    } else if (bits_per_distance == 4) {
      p += (j >> 1);
      distance = *((uint8_t *)p) >> ((j & 1) << 2);
    } else if (bits_per_distance == 8) {
      p += j; 
      distance = *((uint8_t *)p);
    } else if (bits_per_distance == 12) {
      p += j + (j >> 1); 
      distance = *((uint16_t *)p) >> ((j & 1) << 2);
    } else if (bits_per_distance == 16) {
      p += (j << 1); 
      distance = *((uint16_t *)p);
    } else if (bits_per_distance == 32) {
      distance = ((uint32_t *)p)[j];
    }

    return distance & kDistanceMask; 
  }




void PrintAllTags1() const {
    std::cout << "Combined Table contents (Distance, Fingerprint):" << std::endl;

    for (size_t i = 0; i < num_buckets_; ++i) {
        std::cout << "Bucket " << i << ": ";

        for (size_t j = 0; j < kTagsPerBucket; ++j) {
            uint32_t distance = ReadDistance(i, j);
            uint32_t fingerprint = ReadFingerprint(i, j);
            std::cout << std::setw(4) << "(" << distance << "," << fingerprint << ") ";
        }

        std::cout << std::endl;
    }
}


  inline void WriteFingerprint(const size_t i, const size_t j, const uint32_t fp) {

  
    char *p = buckets_fp_[i].bits_; 
    uint32_t fingerprint = fp & kFpMask; 

    if (bits_per_fp == 2) {
      *((uint8_t *)p) &= ~(0x3 << (2 * j)); 
      *((uint8_t *)p) |= fingerprint << (2 * j); 
    } else if (bits_per_fp == 4) {
      p += (j >> 1); 
      if ((j & 1) == 0) {
        *((uint8_t *)p) &= 0xF0; 
        *((uint8_t *)p) |= fingerprint; 
      } else {
        *((uint8_t *)p) &= 0x0F; 
        *((uint8_t *)p) |= (fingerprint << 4); 
      }
    } else if (bits_per_fp == 8) {  
      ((uint8_t *)p)[j] = fingerprint; 
    } else if (bits_per_fp == 12) {
      p += (j + (j >> 1)); 
      if ((j & 1) == 0) {
        ((uint16_t *)p)[0] &= 0xF000; 
        ((uint16_t *)p)[0] |= fingerprint; 
      } else {
        ((uint16_t *)p)[0] &= 0x000F; 
        ((uint16_t *)p)[0] |= (fingerprint << 4); 
      }
    } else if (bits_per_fp == 16) {
      ((uint16_t *)p)[j] = fingerprint; 
    } else if (bits_per_fp == 32) {
      ((uint32_t *)p)[j] = fingerprint; 
    }
    
  }


inline void WriteDistance(const size_t i, const size_t j, const uint32_t distance) {
    char *p = buckets_distance_[i].bits_; 
    uint32_t dist = distance & kDistanceMask; 

    if (bits_per_distance == 2) {
        
        *((uint8_t *)p) &= ~(0x3 << (2 * j)); 
        *((uint8_t *)p) |= dist << (2 * j);   
    } else if (bits_per_distance == 4) {
        
        p += (j >> 1); 
        if ((j & 1) == 0) {
            *((uint8_t *)p) &= 0xF0; 
            *((uint8_t *)p) |= dist; 
        } else {
            *((uint8_t *)p) &= 0x0F; 
            *((uint8_t *)p) |= (dist << 4); 
        }
    } else if (bits_per_distance == 8) {
      
        ((uint8_t *)p)[j] = dist; 
    } else if (bits_per_distance == 16) {
       
        ((uint16_t *)p)[j] = dist; 
    } else if (bits_per_distance == 32) {
        
        ((uint32_t *)p)[j] = dist; 
    } else {
        throw std::invalid_argument("Unsupported bits_per_distance value.");
    }
}


inline bool Find_fp_In_Main_Buckets(const size_t i1, const uint32_t fp, size_t q, size_t p) const {
    for (size_t attempt = 0; attempt < q; attempt++) {
        size_t new_bucket_i = i1 + attempt;
        if (new_bucket_i >= num_buckets_) new_bucket_i -= num_buckets_;

        __m128i target_fp = _mm_set1_epi32(fp);
        __m128i target_dist = _mm_set1_epi32(attempt);

        for (size_t j = 0; j < kTagsPerBucket; j += 4) {
            __m128i fingerprints = _mm_set_epi32(
                ReadFingerprint(new_bucket_i, j + 3),
                ReadFingerprint(new_bucket_i, j + 2),
                ReadFingerprint(new_bucket_i, j + 1),
                ReadFingerprint(new_bucket_i, j)
            );
            __m128i distances = _mm_set_epi32(
                ReadDistance(new_bucket_i, j + 3),
                ReadDistance(new_bucket_i, j + 2),
                ReadDistance(new_bucket_i, j + 1),
                ReadDistance(new_bucket_i, j)
            );

            __m128i fp_cmp = _mm_cmpeq_epi32(fingerprints, target_fp);
            __m128i dist_cmp = _mm_cmpeq_epi32(distances, target_dist);

            __m128i result = _mm_and_si128(fp_cmp, dist_cmp);

            if (_mm_movemask_epi8(result) != 0) {
                return true;
            }
        }
    }
    return false;
}


inline bool Find_fp_In_Auxiliary_Buckets(const size_t i2, const uint32_t fp, size_t q, size_t p) {
    for (size_t attempt = q; attempt < p; attempt++) {
        size_t new_bucket_i = i2 + (attempt - q);
        if (new_bucket_i >= num_buckets_) new_bucket_i -= num_buckets_;

        __m128i target_fp = _mm_set1_epi32(fp);
        __m128i target_dist = _mm_set1_epi32(attempt);

        for (size_t j = 0; j < kTagsPerBucket; j += 4) {
            __m128i fingerprints = _mm_set_epi32(
                ReadFingerprint(new_bucket_i, j + 3),
                ReadFingerprint(new_bucket_i, j + 2),
                ReadFingerprint(new_bucket_i, j + 1),
                ReadFingerprint(new_bucket_i, j)
            );
            __m128i distances = _mm_set_epi32(
                ReadDistance(new_bucket_i, j + 3),
                ReadDistance(new_bucket_i, j + 2),
                ReadDistance(new_bucket_i, j + 1),
                ReadDistance(new_bucket_i, j)
            );

            __m128i fp_cmp = _mm_cmpeq_epi32(fingerprints, target_fp);

            __m128i dist_cmp = _mm_cmpeq_epi32(distances, target_dist);

            __m128i result = _mm_and_si128(fp_cmp, dist_cmp);

            if (_mm_movemask_epi8(result) != 0) {
                return true;
            }
        }
    }
    return false;
}


inline bool DeleteTagFrom_MainBuckets(const size_t i1, const uint32_t fp, size_t q, size_t p) {
    for (size_t attempt = 0; attempt < q; attempt++) {
        size_t new_bucket_i = i1 + attempt;
        if (new_bucket_i >= num_buckets_) new_bucket_i -= num_buckets_;

        for (size_t j = 0; j < kTagsPerBucket; j++) {

            uint32_t current_fp = ReadFingerprint(new_bucket_i, j);
            size_t current_dist = ReadDistance(new_bucket_i, j);

            if (current_fp == fp && current_dist == attempt) {

                WriteFingerprint(new_bucket_i, j, 0);
                WriteDistance(new_bucket_i, j, 0);
                return true; 
            }
        }
    }
    return false; 
}


inline bool DeleteTagFrom_AuxiliaryBuckets(const size_t i2, const uint32_t fp, size_t q, size_t p) {
    for (size_t attempt = q; attempt < p; attempt++) {
        size_t new_bucket_i = i2 + (attempt - q);
        if (new_bucket_i >= num_buckets_) new_bucket_i -= num_buckets_;

        for (size_t j = 0; j < kTagsPerBucket; j++) {

            uint32_t current_fp = ReadFingerprint(new_bucket_i, j);
            size_t current_dist = ReadDistance(new_bucket_i, j);

            if (current_fp == fp && current_dist == attempt) {

                WriteFingerprint(new_bucket_i, j, 0);
                WriteDistance(new_bucket_i, j, 0);
   
                return true; 
            }
        }
    }
    return false; 
}




inline size_t NumTagsInBucket(const size_t i) const {
    size_t num = 0;
    for (size_t j = 0; j < kTagsPerBucket; j++) {
        if (ReadFingerprint(i, j) != 0) {
            num++; 
        }
    }
    return num;
}

};

} // namespace HCF
#endif  //HCF
