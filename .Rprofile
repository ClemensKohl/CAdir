source("renv/activate.R")
options(renv.download.override = utils::download.file)


# General ----------------------------------------
if (require("colorout", quietly = TRUE)) {

    colorout::setOutputColors(

                              index    = "\x1b[38;2;115;121;148m", #'\x1b[38;2;76;86;106m',
                              normal   = "\x1b[38;2;140;170;238m", #"\x1b[38;2;198;208;245m", #'\x1b[38;2;216;222;233m',

                              number   = "\x1b[38;2;242;213;207m",  #'\x1b[38;2;236;239;244m',
                              negnum   = "\x1b[38;2;244;184;228m", #'\x1b[38;2;180;142;173m',
                              zero     = "\x1b[38;2;238;190;190m", #"\x1b[38;2;186;187;241m", #'\x1b[38;2;136;192;208m',
                              zero.limit = 0.01,
                              infinite = "\x1b[38;2;244;184;228m", #'\x1b[38;2;236;239;244m',

                              string   = "\x1b[38;2;166;209;137m", #'\x1b[38;2;235;203;139m',
                              date     = "\x1b[38;2;229;200;144m", #'\x1b[38;2;236;239;244m',
                              const    = "\x1b[38;2;229;200;144m", #'\x1b[38;2;136;192;208m',

                              true     = "\x1b[38;2;129;200;190m", #'\x1b[38;2;163;190;140m',
                              false    = "\x1b[38;2;244;184;228m", #'\x1b[38;2;191;97;106m',

                              warn     = "\x1b[38;2;239;159;118m", #'\x1b[38;2;235;203;139m',
                              stderror = "\x1b[38;2;186;187;241m", #"\x1b[38;2;231;130;132m", #'\x1b[38;2;191;97;106m',
                              error = "\x1b[38;2;231;130;132m", #'\x1b[38;2;191;97;106m',

                              verbose  = FALSE
    )

    # Class
    colorout::addPattern(' num ',        "\x1b[38;2;129;200;190m")
    colorout::addPattern(' int ',        "\x1b[38;2;129;200;190m")
    colorout::addPattern(' chr ',        "\x1b[38;2;129;200;190m")
    colorout::addPattern(' Factor ',     "\x1b[38;2;129;200;190m")
    colorout::addPattern(' Ord.factor ', "\x1b[38;2;129;200;190m")
    colorout::addPattern(' logi ',       "\x1b[38;2;129;200;190m")
    colorout::addPattern('function ',    "\x1b[38;2;129;200;190m")
    colorout::addPattern(' dbl ',        "\x1b[38;2;129;200;190m")
    colorout::addPattern(' lgcl ',       "\x1b[38;2;129;200;190m")
    colorout::addPattern(' cplx ',       "\x1b[38;2;129;200;190m")
    # Misc
    colorout::addPattern('$ ',           "\x1b[38;2;129;200;190m")

    # _ `str`, {mlr3} --------------------------------

    # R6 field name
    colorout::addPattern('* [A-z]*:',                      '\x1b[38;2;235;203;139m')
    colorout::addPattern("* [A-z]* [A-z]*:",               '\x1b[38;2;235;203;139m')
    colorout::addPattern("* [A-z]* [A-z]* [A-z]*:",        '\x1b[38;2;235;203;139m')
    colorout::addPattern("* [A-z]* [A-z]* [A-z]* [A-z]*:", '\x1b[38;2;235;203;139m')
    # So on...

    #############

}
