# Prove the tree with the switches the units were admitted under.
# Requires gnatprove FSF 15.0 on PATH (see README.md).

GNATPROVE ?= gnatprove
JOBS ?= 4

.PHONY: prove clean checksums

prove:
	$(GNATPROVE) -P pde_kernels.gpr --level=2 --report=all -j$(JOBS)

checksums:
	shasum -a 256 -c SHA256SUMS

clean:
	$(GNATPROVE) -P pde_kernels.gpr --clean
	rm -rf obj
