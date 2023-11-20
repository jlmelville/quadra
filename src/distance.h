// BSD 2-Clause License
//
// Copyright 2022 James Melville
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// OF SUCH DAMAGE.

#ifndef QUADRA_DISTANCE_H
#define QUADRA_DISTANCE_H

#include <functional>
#include <string>

#include "tdoann/distance.h"

using It = typename std::vector<double>::const_iterator;

inline std::function<double(It, It, It)>
create_dfun(const std::string &metric) {
  if (metric == "euclidean") {
    return tdoann::euclidean<double, It>;
  } else if (metric == "sqeuclidean") {
    return tdoann::squared_euclidean<double, It>;
  } else if (metric == "cosine") {
    return tdoann::cosine<double, It>;
  } else if (metric == "hamming") {
    return tdoann::hamming<double, It>;
  } else if (metric == "manhattan") {
    return tdoann::manhattan<double, It>;
  } else if (metric == "correlation") {
    return tdoann::correlation<double, It>;
  } else {
    return tdoann::euclidean<double, It>;
  }
}

#endif // QUADRA_DISTANCE_H
