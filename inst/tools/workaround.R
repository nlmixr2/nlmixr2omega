## This is only for nlmixr2omega

.in <- suppressWarnings(readLines("src/RcppExports.cpp"))
.w <- which(regexpr("Rcpp[.]h", .in) != -1)[1]
if (regexpr("nlmixr2omegaArma[.]h", .in[.w+1]) == -1) {
  .in <- c(.in[seq(1,.w)],
           '#include "../inst/include/nlmixr2omegaArma.h"',
           .in[seq(.w+1, length(.in))])
  file.out <- file("src/RcppExports.cpp", "wb")
  writeLines(.in, file.out)
  close(file.out)
}


.in <- suppressWarnings(readLines("src/Makevars.in"))
.in <- gsub("@ARMA@", file.path(find.package("RcppArmadillo"),"include"), .in)
#.in <- gsub("@BH@", file.path(find.package("BH"),"include"), .in)
.in <- gsub("@RCPP@", file.path(find.package("Rcpp"),"include"), .in)
#.in <- gsub("@EG@", file.path(find.package("RcppEigen"),"include"), .in)

.badStan <- ""
.in <- gsub("@SH@", gsub("-I", "-@ISYSTEM@",
                         paste(## capture.output(StanHeaders:::CxxFlags()),
                               ## capture.output(RcppParallel:::CxxFlags()),
                               paste0("-@ISYSTEM@'", system.file('include', 'src', package = 'StanHeaders', mustWork = TRUE), "'"),
                               .badStan)),
            .in)

.in <- gsub("@SL@", "", ##paste(capture.output(StanHeaders:::LdFlags()), capture.output(RcppParallel:::RcppParallelLibs())),
            .in)


if (.Platform$OS.type == "windows" && !file.exists("src/Makevars.win")) {
  .in <- gsub("@CXX14STD@", "-std=c++1y", .in)
  file.out <- file("src/Makevars.win", "wb")
  writeLines(gsub("@ISYSTEM@", "I", .in),
             file.out)
  close(file.out)
} else {
  .in <- gsub("@CXX14STD@", "-std=gnu++14", .in)
  file.out <- file("src/Makevars", "wb")
  writeLines(gsub("@ISYSTEM@", "isystem", .in),
             file.out)
  close(file.out)
}

unlink("R/nlmixr2omega_md5.R")

cpp <- list.files("src", pattern = ".(c|h|cpp|f)$")
include <- list.files("inst/include")
#Rfiles <- list.files("R/", pattern = ".R")

cmd <- file.path(R.home("bin"), "R")
args <- c("CMD", "config")

md5 <- digest::digest(c(lapply(c(paste0("src/", cpp),
                                 paste0("inst/include/", include)#,
                                 #paste0("R/", Rfiles)
                                 ), digest::digest, file = TRUE),
                        ## vapply(c("BLAS_LIBS", "CC",  "CFLAGS", "CPICFLAGS",
                        ##          "CXX", "CXXFLAGS", "CXXPICFLAGS",
                        ##          "CXX11", "CXX11STD", "CXX11FLAGS", "CXX11PICFLAGS",
                        ##          "CXX14", "CXX14STD", "CXX14FLAGS", "CXX14PICFLAGS",
                        ##          "CXX17", "CXX17STD", "CXX17FLAGS", "CXX17PICFLAGS",
                        ##          "CXX20", "CXX20STD", "CXX20FLAGS", "CXX20PICFLAGS",
                        ##          "FC", "FFLAGS", "FCFLAGS",  "FPICFLAGS"),
                        ##        function(cfg) {
                        ##          rawToChar(sys::exec_internal(cmd, c(args, cfg))$stdout)
                        ##        }, character(1)
                        ##       ),
                        ""
                        ))
unlink("R/nlmixr2omega_md5.R")
md5file <- file("R/nlmixr2omega_md5.R", "wb")
writeLines(sprintf("nlmixr2omega.md5 <- \"%s\"\n", md5), md5file)
close(md5file)
