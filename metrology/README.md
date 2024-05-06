
## LVM Fibermap

This readme describes the LVM metrology files, how to update them, and how

The main LVM fiber mapping is in the `lvm_fiducial_fibermap.yaml` file.

It contains two sections:
- **schema**: a description of the columns in the fibermap
- **fibers**: a list of fiber information, with a fiber per row

### Swapping Fibers

If fibers are discovered to be swapped, the rows in the fibermap file need updating. To properly swap the fibers, update the following columns with the correct information:

- **fiberid**
- **finblock**
- the slitblock info in **fmap**

Then **swap the rows** in the fibermap file.  This is to preserve the fibermap order so it properly matches the reduced data, without requiring a manual sort.

**Warning:** The `ypix_x` (where x is either 'b', 'r' or 'z') column should never be swapped, as this is position, not fiber, dependent.  When moving rows around, remember to keep the `ypix_x` values in their
original positions.

**Example**

**Old**
Fibers 431/432 are incorrect. The original S1-400 science fiber and the B1-13 SKY fiber, should be swapped.
```
  - [430, 1, B12, 34, science, S1, 399, Sci, 4.001, 0.0, 15.0, S1-399, S1B12-34, 399, 'Sci1-399:S1B12-34', 1476, 1476, 1469, 0]
  - [431, 1, B12, 35, science, S1, 400, Sci, 4.001, -0.33, 15.0, S1-400, S1B12-35, 400, 'Sci1-400:S1B12-35', 1471, 1471, 1463, 0]
  - [432, 1, B12, 36, SKY, B1, 13, SkyE, 0.572, 0.99, 5.0, B1-12, S1B12-36, 12, 'SkyE-12:S1B12-36', 1465, 1465, 1458, 0]
  - [433, 1, B13, 1, SKY, A1, 14, SkyW, 0.857, 0.825, 5.0, A1-13, S1B13-1, 13, 'SkyW-13:S1B13-1', 1448, 1447, 1440, 0]
```
**New**
The fiberid, finblock, and fmap columns have been updated. The rows have been reordered to reflect the change.  The ypix_x column remains the same.
```
  - [430, 1, B12, 34, science, S1, 399, Sci, 4.001, 0.0, 15.0, S1-399, S1B12-34, 399, 'Sci1-399:S1B12-34', 1476, 1476, 1469, 0]
  - [431, 1, B12, 35, SKY, B1, 13, SkyE, 0.572, 0.99, 5.0, B1-12, S1B12-36, 12, 'SkyE-12:S1B12-35', 1471, 1471, 1463, 0]
  - [432, 1, B12, 36, science, S1, 400, Sci, 4.001, -0.33, 15.0, S1-400, S1B12-35, 400, 'Sci1-400:S1B12-36', 1465, 1465, 1458, 0]
  - [433, 1, B13, 1, SKY, A1, 14, SkyW, 0.857, 0.825, 5.0, A1-13, S1B13-1, 13, 'SkyW-13:S1B13-1', 1448, 1447, 1440, 0]
```

### Updating Fiber Status

The last column is the fiber status.  It is an integer value with the following meanings:

- 0 - ok: good fiber
- 1 - dead: dead fibers
- 2 - low: fibers with visually low throughput
- 3 - repair: repaired fibers
- 4 - short: fiber short at splice end; visible at IFU end

To update the fiber status, change the integer value in the `fibstatus` column.

