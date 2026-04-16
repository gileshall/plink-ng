#ifndef __PLINK2_HAP_IBD_H__
#define __PLINK2_HAP_IBD_H__

// This file is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "include/pgenlib_misc.h"
#include "include/pgenlib_read.h"
#include "include/plink2_base.h"
#include "plink2_cmdline.h"
#include "plink2_common.h"

#ifdef __cplusplus
namespace plink2 {
#endif

ENUM_U31_DEF_START()
  kHapIbdMethodPbwt,
  kHapIbdMethodPbwtAdaptive,
  kHapIbdMethodPbwtHmm
ENUM_U31_DEF_END(HapIbdMethod);

FLAGSET_DEF_START()
  kfHapIbd0,
  kfHapIbdExtendDiploid = (1 << 0),
  kfHapIbdWhichIbd1 = (1 << 1),
  kfHapIbdWhichIbd2 = (1 << 2),
  kfHapIbdWhichBoth = (kfHapIbdWhichIbd1 | kfHapIbdWhichIbd2),
  kfHapIbdOutFmtSegments = (1 << 3),
  kfHapIbdOutFmtSummary = (1 << 4),
  kfHapIbdOutFmtMatrix = (1 << 5),
  kfHapIbdOutFmtSparse = (1 << 6),
  kfHapIbdOutFmtMask = (kfHapIbdOutFmtSegments | kfHapIbdOutFmtSummary | kfHapIbdOutFmtMatrix | kfHapIbdOutFmtSparse),
  kfHapIbdDiag = (1 << 7)
FLAGSET_DEF_END(HapIbdFlags);

typedef struct HapIbdInfoStruct {
  NONCOPYABLE(HapIbdInfoStruct);
  HapIbdFlags flags;
  HapIbdMethod method;
  uint32_t min_l_bp;       // --hap-ibd-min-l (default 1000000)
  uint32_t min_snp;        // --hap-ibd-min-snp (default 100)
  uint32_t max_gap;        // --hap-ibd-max-gap in bp (default 500000)
  uint32_t max_err;        // --hap-ibd-max-err (default 0 for pbwt)
  uint32_t seed_len;       // --hap-ibd-seed-len (default 50)
  uint32_t trim_bp;        // --hap-ibd-trim (default 0)
  uint32_t seg_buf_size;   // --hap-ibd-seg-buf in MB (default 256)
  uint32_t thread_ct;      // --hap-ibd-threads (0 = use global)
  double min_cm;           // --hap-ibd-min-cm (default -1 = unset)
  double err_rate;         // --hap-ibd-err-rate (default 0.002)
  double maf_cap;          // --hap-ibd-maf-cap (default 1.0 = no cap)
  double min_kin;          // --hap-ibd-min-kin (default 0.01)
  char* cm_map_fname;      // --hap-ibd-cm genetic map file
} HapIbdInfo;

void InitHapIbd(HapIbdInfo* hap_ibd_info_ptr);

void CleanupHapIbd(HapIbdInfo* hap_ibd_info_ptr);

PglErr CalcHapIbd(const uintptr_t* sample_include, const SampleIdInfo* siip, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const double* variant_cms, const double* allele_freqs, const HapIbdInfo* hap_ibd_info_ptr, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_thread_ct, PgenFileInfo* pgfip, PgenReader* simple_pgrp, char* outname, char* outname_end);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_HAP_IBD_H__
