# esreg 1.0.0 -- the copy submitted to the SSC archive

The same files as esreg 1.0.0 on GitHub (tag v1.0.0), with the technical note as
an ancillary file (esreg_technical_note.pdf). It is here to test the SSC
package before the archive publishes it. As on SSC, the package file lists
every file with an `f` line: `net install` installs the programs, the help
and the dialog only; the example data and the technical note are copied by
`net get` to the current folder, never to the system directories.

```stata
net install esreg, from("https://raw.githubusercontent.com/aabbdd12/esreg/main/ssc/1.0.0") replace
net get esreg, from("https://raw.githubusercontent.com/aabbdd12/esreg/main/ssc/1.0.0") replace
```

The examples of the help read the data from the current folder, else from the
SSC archive, else from GitHub. To go back to the GitHub version:

```stata
net install esreg, from("https://raw.githubusercontent.com/aabbdd12/esreg/v1.0.0") replace
```
