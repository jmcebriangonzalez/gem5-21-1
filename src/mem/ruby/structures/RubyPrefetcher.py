# Copyright (c) 2020 ARM Limited
# All rights reserved.
#
# The license below extends only to copyright in the software and shall
# not be construed as granting a license to any other intellectual
# property including but not limited to intellectual property relating
# to a hardware implementation of the functionality of the software
# licensed hereunder.  You may use the software subject to the license
# terms below provided that you ensure that this notice is replicated
# unmodified and in its entirety in all distributions of the software,
# modified or unmodified, in source code or in binary form.
#
# Copyright (c) 2012 Mark D. Hill and David A. Wood
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met: redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer;
# redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution;
# neither the name of the copyright holders nor the names of its
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from m5.SimObject import SimObject
from m5.params import *
from m5.proxy import *

from m5.objects.System import System

class RubyPrefetcher(SimObject):
    type = 'RubyPrefetcher'
    cxx_class = 'gem5::ruby::RubyPrefetcher'
    cxx_header = "mem/ruby/structures/RubyPrefetcher.hh"

    log_num_sets_in_recent_access_tag_array_l1 = Param.UInt64(0, "")
    num_ways_in_recent_access_tag_array_l1 = Param.UInt64(40, "")

    ip_table_tag_mask = Param.UInt64(0x3fff, "Prefetch degree 4")
    log_num_sets_in_ip_table_l1 = Param.UInt64(7, "")

    num_ways_in_ip_table_l1 = Param.UInt64(15, "")
    ip_delta_table_tag_mask = Param.UInt64(0xffff, "Prefetch degree 4")
    log_num_sets_in_ip_delta_table_l1 = Param.UInt64(8, "")
    num_ways_in_ip_delta_table_l1 = Param.UInt64(8, "")

    saturating_counter_max_l1 = Param.UInt64(3, "")
    base_prefetch_degree_l1 = Param.UInt64(4, "")
    num_entries_in_nl_buffer_l1 = Param.UInt64(64, "")
    nl_threshold_numer_l1 = Param.UInt64(1, "")
    nl_threshold_denom_l1 = Param.UInt64(4, "")

    pointer_last = Param.Bool(True, "")
    pointer_non_last = Param.Bool(False, "")
    stride_conf_max = Param.UInt64(3, "")
    stride_conf_threshold = Param.UInt64(3, "")

    partial_ip_mask = Param.UInt64(0x7f, "")

    num_strides_in_long_hist_ip_table = Param.UInt64(20, "")
    long_hist_ip_table_tag_mask = Param.UInt64(0x1fffff, "")
    num_entries_in_long_hist_ip_table = Param.UInt64(32, "")
    long_hist_match_length = Param.UInt64(1, "")

    num_ip_table_l1_entries = Param.UInt64(1024, "ip table entries")
    num_ghb_entries = Param.UInt64(16, "entries in the ghb")
    num_ip_index_bits = Param.UInt64(10, "bits to index into the ip table")
    num_ip_tag_bits = Param.UInt64(6, "tag bits per ip table entry")
    s_type = Param.UInt64(1, "stream")
    cs_type = Param.UInt64(2, "constant stride")
    cplx_type = Param.UInt64(3, "complex stride")
    nl_type = Param.UInt64(4, "next line")

    sys = Param.System(Parent.any, "System this prefetcher belongs to")
class Prefetcher(RubyPrefetcher):
    """DEPRECATED"""
    pass
