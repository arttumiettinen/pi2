.. _readmeta:

readmeta
********


**Syntax:** :code:`readmeta(image, filename)`

Reads image metadata from disk.

This command cannot be used in the distributed processing mode. If you need it, please contact the authors.

Arguments
---------

image [input & output]
~~~~~~~~~~~~~~~~~~~~~~

**Data type:** uint8 image, uint16 image, uint32 image, uint64 image, int8 image, int16 image, int32 image, int64 image, float32 image, complex32 image

Image to process.

filename [input]
~~~~~~~~~~~~~~~~

**Data type:** string

Name and path of file to read from.

See also
--------

:ref:`setmeta`, :ref:`getmeta`, :ref:`writemeta`, :ref:`readmeta`, :ref:`clearmeta`, :ref:`listmeta`
