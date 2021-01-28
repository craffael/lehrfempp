/**
 * @file
 * @brief Implementation of LehrFemInfo class.
 * Note that some information in this file substituted by cmake during build
 * using cmake's `configure_file` command.
 * @author Raffael Casagrande
 * @date   2020-10-23 09:08:09
 * @copyright MIT License
 */

#include <lf/base/lehrfem_info.h>

#include <regex>
#include <sstream>

namespace /*anonymous*/ {}

namespace lf::base {
std::string LehrFemInfo::getVersionSha() { return "@version_sha@"; }

std::string LehrFemInfo::getVersionDateTime() { return "@version_datetime@"; }

std::string LehrFemInfo::getVersionTag() { return "@version_tag@"; }

void LehrFemInfo::PrintInfo(std::ostream &stream) {
  std::string version;
  if (!getVersionTag().empty()) {
    version = getVersionTag();
  } else if (!getVersionSha().empty()) {
    version = getVersionSha();
  } else {
    version = "unknown";
  }
  version.insert(version.size(), 56 - version.size(), ' ');
  version += "|\n";

  std::string version_datetime;
  if (!getVersionDateTime().empty()) {
    version_datetime = "(" + getVersionDateTime() + ")";
  }

  version_datetime.insert(version_datetime.size(), 56 - version_datetime.size(),
                          ' ');
  version_datetime += "|\n";

  // clang-format off
  stream << "\n";
  stream << R"foo(+------------------------------------------------------------------+)foo" << "\n";
  stream << R"foo(|        __         __         ______________  ___                 |)foo" << "\n";
  stream << R"foo(|       / /   ___  / /_  _____/ ____/ ____/  |/  / __    __        |)foo" << "\n";
  stream << R"foo(|      / /   / _ \/ __ \/ ___/ /_  / __/ / /|_/ /_/ /___/ /_       |)foo" << "\n";
  stream << R"foo(|     / /___/  __/ / / / /  / __/ / /___/ /  / /_  __/_  __/       |)foo" << "\n";
  stream << R"foo(|    /_____/\___/_/ /_/_/  /_/   /_____/_/  /_/ /_/   /_/          |)foo" << "\n";
  stream << R"foo(|                                                                  |)foo" << "\n";
  stream << R"foo(|                                                                  |)foo" << "\n";
  stream << R"foo(| Maintained by Raffael Casagrande and Ralf Hiptmair with          |)foo" << "\n";
  stream << R"foo(| contributions from many other people.                            |)foo" << "\n";
  stream << R"foo(|                                                                  |)foo" << "\n";
  stream << R"foo(| Version: )foo" << version;           
  stream << R"foo(|          )foo" << version_datetime;
  stream << R"foo(|  Source: https://github.com/craffael/lehrfempp                   |)foo" << "\n";
  stream << R"foo(| License: MIT License                                             |)foo" << "\n";
  stream << R"foo(+------------------------------------------------------------------+)foo" << "\n";
  stream << "\n";
  // clang-format on
}

void LehrFemInfo::PrintLicense(std::ostream &stream) {
  stream << R"foo(MIT License

Copyright (c) 2018 Raffael Casagrande, Ralf Hiptmair

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

)foo";
}

void LehrFemInfo::Print3rdPartyLicenses(std::ostream &stream) {
  stream << R"foo(
LehrFEM++ makes use of the following 3rd Party software with their respective License Terms:

Hunter Package Manager (https://github.com/cpp-pm/hunter)
  BSD 2-Clause "Simplified License"
---------------------------------------------------------
Copyright (c) 2013-2018, Ruslan Baratov
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Boost (https://boost.org)
---------------------------------------------------------
Boost Software License - Version 1.0 - August 17th, 2003

Permission is hereby granted, free of charge, to any person or organization
obtaining a copy of the software and accompanying documentation covered by
this license (the "Software") to use, reproduce, display, distribute,
execute, and transmit the Software, and to prepare derivative works of the
Software, and to permit third-parties to whom the Software is furnished to
do so, all subject to the following:

The copyright notices in the Software and this entire statement, including
the above license grant, this restriction and the following disclaimer,
must be included in all copies of the Software, in whole or in part, and
all derivative works of the Software, unless such copies or derivative
works are solely in the form of machine-executable object code generated by
a source language processor.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.

Eigen (http://eigen.tuxfamily.org/)
---------------------------------------------------------
Eigen uses the Mozilla Public License Version 2.0
(http://mozilla.org/MPL/2.0/)

Eigens source code is available at https://gitlab.com/libeigen/eigen

fmt (https://github.com/fmtlib/fmt)
---------------------------------------------------------
Copyright (c) 2012 - present, Victor Zverovich

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

--- Optional exception to the license ---

As an exception, if, as a result of your compiling your source code, portions
of this Software are embedded into a machine-executable object form of such
source code, you may redistribute such embedded portions in such object form
without including the above copyright and permission notices.


spdlog (https://github.com/gabime/spdlog)
---------------------------------------------------------
The MIT License (MIT)

Copyright (c) 2016 Gabi Melman.                                       

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

Google Test (https://github.com/google/googletest)
---------------------------------------------------------
Copyright 2008, Google Inc.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above
copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the
distribution.
    * Neither the name of Google Inc. nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

)foo";
}

}  // namespace lf::base
