# Quick Reference - Logging {#quick_reference_logging}

[TOC]

## Overview

LehrFem relies on the [spdlog](https://github.com/gabime/spdlog) library for logging.

```cpp
#include "lf/base/base.h"

int main()
{
    // Create a logger, @ref lf::base::InitLogger creates a logger with sensible defaults
    auto logger = lf::base::InitLogger("lf::mesh::hybrid2d::Mesh::Logger");
    logger->info("Hello, World!");
    // Will output: [Time & Date] [info] Hello, World!
    logger->warn("This is a warning!");
    // Will output: [Time & Date] [warn] This is a warning!
    logger->error("This is an error!");
    // Will output: [Time & Date] [error] This is an error!

    return 0;
}
```

The `InitLogger` function creates a logger with the given name.

# Logging Levels

Setting the logging level determines which messages are logged. Log messages with a level lower than the set level are ignored.

The log levels in order of increasing severity are:

-   trace `spdlog::level::trace`
-   debug `spdlog::level::debug`
-   info `spdlog::level::info`
-   warn `spdlog::level::warn`
-   err `spdlog::level::err`
-   critical `spdlog::level::critical`

By default, the logging level is set to `spdlog::level::info`, which means that only messages with level above or equal to `info` are logged (`info`, `warn`, `err`, `critical`). The following example shows how to set the logging level to `trace`:

```cpp
logger->set_level(spdlog::level::trace);
```

# Custom Logging Configuration

Spdlog can be customized further (color, files, etc.) and detailed information and more examples can be found in the [spdlog Github repo](https://github.com/gabime/spdlog)

# Logging in LehrFEM++

Many classes in LehrFEM++ use the logging functionality. For example, the `Mesh` class logs information about the mesh when it is created. The following example shows how to create a mesh and log information about it:

```cpp

```

## Eigen Matrix Logging

Eigen matrices are pretty-printed when logged. The following example shows how to log an Eigen matrix:

```cpp
Eigen::MatrixXd A(2, 2);
A << 1, 2, 3, 4;

// output the matrix using the spdlog logger
logger->info("Matrix A = \n{}", A);
```

The output will look like this:

```console
[Time & Date] info] Matrix A =
                            [1, 2]
                            [3, 4]

```

# Logging in LehrFEM++ Examples

TODO: Add LehrFEM++ example that uses logging
