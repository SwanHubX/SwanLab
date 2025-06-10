<div align="center">
    <picture>
      <source media="(prefers-color-scheme: dark)" srcset="../readme_files/swanlab-logo-single-dark.svg">
      <source media="(prefers-color-scheme: light)" srcset="../readme_files/swanlab-logo-single.svg">
      <img alt="SwanLab" src="../readme_files/swanlab-logo-single.svg" width="70" height="70">
    </picture>
</div>

# swanlab-core: A new backend for SwanLab SDK

To address the challenges posed by Python's asynchronous programming, we introduced swanlab-core, a brand-new, Go-based backend for metric uploading.


## IDE & Development Setup

In this part, we will guide you through the process of setting up your IDE and development environment for SwanLab SDK.
We will use:
1. [GoLand](https://www.jetbrains.com/go/) as the IDE.
2. [Go](https://go.dev/) as the programming language.
3. [golangci-lint](https://golangci-lint.run/) as the linter.

## IDE

We request using GoLand as the IDE for swanlab-core development. 
It provides a rich set of features for Go development, including code completion, debugging, and integration with version control systems.

There is the code template for Go file in GoLand:

```text
package ${GO_PACKAGE_NAME}

// @Title        ${FILE_NAME}
// @Description  
// @Create       cunyue ${DATE} ${TIME}
```

## Linter

We use [golangci-lint](https://golangci-lint.run/) as the linter for swanlab-core development.
Since SwanLab is a Python project, some developers may not be very familiar with the Go language. 
Therefore, golangci-lint is not mandatory and will only trigger this CI when changes are made to the core part.

If you plan to develop swanlab-core, configure the golangci-lint integration in GoLand, which can be found under Go -> Linter settings.

## Env

TODO

