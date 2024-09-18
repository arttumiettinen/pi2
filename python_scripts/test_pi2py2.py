
import pytest
import os
import shutil
import sys

from pi2py2 import *
import numpy as np
import zarrita


# Define convenience functions that return input and output file names.
# This is just to avoid copying the paths to all the examples in case they change.
def input_file(filename='t1-head_256x256x129.raw'):
    return '../test_input_data/' + filename

def input_file_bin():
    return '../test_input_data/t1-head_bin_256x256x129.raw'

def output_file(name):
    return '../test_output_data/pi2py2/' + name

class Test_zarr:

    pi2 = Pi2(library_path='../x64/Release no OpenCL')

    def output_file(self, name):
        return "testoutput/" + name

    fill_value = 42
    full_arr = np.arange(10 * 10 * 10, dtype=np.float32).reshape(10, 10, 10)
    w = 2
    h = 3
    d = 5
    arr = full_arr[:w, :h, :d]
    arr0 = np.zeros_like(arr)
    arr42 = np.zeros_like(arr)
    arr42[:, :, 1] = fill_value

    def pi2_write(self, chunk_shape=None, codecs=None, data=None, name="test.zarr"):
        pi2 = self.pi2

        if data is None:
            data = self.arr
        if chunk_shape is None:
            chunk_shape = list(data.shape)
        shutil.rmtree(output_file(name), ignore_errors=True)
        write_img = pi2.newimage(pi2py2.ImageDataType.UInt32, list(data.shape))
        write_img.set_data(data.transpose(1, 0, 2))
        if codecs is not None:
            pi2.writezarr(write_img, output_file(name), chunk_shape, codecs)
        else:
            pi2.writezarr(write_img, output_file(name), chunk_shape)

    def zarrita_write(self, chunk_shape=None, codecs=None, data=None, separator=None, name="zarrita.zarr"):
        if data is None:
            data = self.arr
        if chunk_shape is None:
            chunk_shape = data.shape
        if separator is None:
            separator = "/"
        if codecs is None:
            codecs = [
                zarrita.codecs.bytes_codec("little"),
            ]

        shutil.rmtree(output_file(name), ignore_errors=True)
        store = zarrita.LocalStore(output_file(name))
        a = zarrita.Array.create(
            store,
            shape=data.shape,
            dtype='float32',
            chunk_shape=tuple(chunk_shape),
            fill_value=self.fill_value,
            chunk_key_encoding=("default", separator),
            codecs=codecs,
        )
        a[:] = data

    def pi2_read(self, name):
        pi2 = self.pi2

        read_img = pi2.read(output_file(name))
        read_arr = read_img.get_data()
        return read_arr.transpose(1, 0, 2)

    def zarrita_read(self, name):
        store = zarrita.LocalStore(output_file(name))
        a = zarrita.Array.open(store)
        read_arr = a[:]
        return read_arr

    @pytest.mark.parametrize("chunk_shape", [[1, 1, 1], [1, 1, d], [1, h, d], [w, h, d]])
    def test_pi2_to_pi2(self, chunk_shape):
        self.pi2_write(chunk_shape)
        read_arr = self.pi2_read("test.zarr")
        assert np.array_equal(self.arr, read_arr), "read_arr:\n " + str(read_arr) + " \n\narr:\n " + str(self.arr)

    @pytest.mark.parametrize("chunk_shape", [[1, 1, 1], [1, 1, d], [1, h, d], [w, h, d]])
    def test_zarrita_to_zarrita(self, chunk_shape):
        self.zarrita_write(chunk_shape)
        read_arr = self.zarrita_read("zarrita.zarr")
        assert np.array_equal(self.arr, read_arr), "read_arr:\n " + str(read_arr) + " \n\narr:\n " + str(self.arr)

    @pytest.mark.parametrize("chunk_shape", [[1, 1, 1], [1, 1, d], [1, h, d], [w, h, d]])
    def test_pi2_to_zarrita(self, chunk_shape):
        self.pi2_write(chunk_shape)
        read_arr = self.zarrita_read("test.zarr")
        assert np.array_equal(self.arr, read_arr), "read_arr:\n " + str(read_arr) + " \n\narr:\n " + str(self.arr)

    @pytest.mark.parametrize("chunk_shape", [[1, 1, 1], [1, 1, d], [1, h, d], [w, h, d]])
    def test_zarrita_to_pi2(self, chunk_shape):
        self.zarrita_write(chunk_shape)
        read_arr = self.pi2_read("zarrita.zarr")
        # write_img = pi2.newimage(pi2py2.ImageDataType.FLOAT32, w, h, d)
        # write_img.set_data(read_arr)
        # pi2.writezarr(write_img, output_file("test.zarr"), chunk_shape)
        assert np.array_equal(self.arr, read_arr), "read_arr:\n " + str(read_arr) + " \n\narr:\n " + str(self.arr)

    @pytest.mark.parametrize("chunk_shape", [[1, 1, 1], [1, 1, d], [1, h, d], [w, h, d]])
    @pytest.mark.parametrize("order", [[0, 1, 2], [1, 0, 2], [0, 2, 1], [1, 2, 0], [2, 0, 1], [2, 1, 0], ])
    def test_read_transpose(self, order, chunk_shape):
        self.zarrita_write(codecs=[zarrita.codecs.transpose_codec(order), zarrita.codecs.bytes_codec("little")],
                      chunk_shape=chunk_shape)
        read_arr = self.pi2_read("zarrita.zarr")
        assert np.array_equal(self.arr, read_arr), "read_arr:\n " + str(read_arr) + " \n\narr:\n " + str(self.arr)

    @pytest.mark.parametrize("cname", ["lz4", "lz4hc", "blosclz", "zstd", "snappy", "zlib"])
    @pytest.mark.parametrize("chunk_shape", [[1, 1, 1], [1, 1, d], [1, h, d], [w, h, d]])
    @pytest.mark.parametrize("typesize", [1, 2, 10])
    @pytest.mark.parametrize("data", [arr, np.zeros_like(arr)])
    # todo: blocksize
    def test_read_blosc(self, cname, chunk_shape, typesize, data):
        self.zarrita_write(codecs=[zarrita.codecs.bytes_codec("little"), zarrita.codecs.blosc_codec(typesize)],
                      chunk_shape=chunk_shape, data=data)
        read_arr = self.pi2_read("zarrita.zarr")
        assert np.array_equal(data, read_arr), "read_arr:\n " + str(read_arr) + " \n\narr:\n " + str(data)

    @pytest.mark.parametrize("codecs", [None, "[{\"configuration\": {\"endian\": \"little\"},\"name\": \"bytes\"}]",
                                        '[{"configuration": {"endian": "little"},"name": "bytes"}]'])
    def test_api_parse_codecs(self, codecs):
        self.pi2_write(codecs=codecs)
        read_arr = self.zarrita_read("test.zarr")
        assert np.array_equal(self.arr, read_arr), "read_arr:\n " + str(read_arr) + " \n\narr:\n " + str(self.arr)

    @pytest.mark.parametrize("chunk_shape", [[1, 1, 1], [1, 1, d], [1, h, d], [w, h, d]])
    @pytest.mark.parametrize("order", [[0, 1, 2], [1, 0, 2], [0, 2, 1], [1, 2, 0], [2, 0, 1], [2, 1, 0], ])
    def test_write_transpose(self, order, chunk_shape):
        codecs = '[{"configuration": {"order": ' + str(
            order) + '},"name": "transpose"}, {"configuration": {"endian": "little"},"name": "bytes"}]'
        self.pi2_write(codecs=codecs, chunk_shape=chunk_shape)
        read_arr = self.zarrita_read("test.zarr")
        assert np.array_equal(self.arr, read_arr), "read_arr:\n " + str(read_arr) + " \n\narr:\n " + str(self.arr)

    @pytest.mark.parametrize("chunk_shape", [[5, 5, 2], [10, 1, 10], [1, 10, 10], [10, 10, 10]])
    @pytest.mark.parametrize("cname", ["lz4", "lz4hc", "blosclz", "zstd", "zlib"])
    @pytest.mark.parametrize("clevel", [1, 2, 4])
    @pytest.mark.parametrize("shuffle", ["shuffle"])
    @pytest.mark.parametrize("typesize", [4])
    @pytest.mark.parametrize("blocksize", [0])
    @pytest.mark.parametrize("data", [full_arr, np.zeros_like(full_arr)])
    # @pytest.mark.parametrize("cname", ["lz4", "lz4hc", "blosclz", "zstd", "snappy", "zlib"])
    def test_write_blosc(self, chunk_shape, cname, clevel, shuffle, typesize, blocksize, data):
        codecs = '[{"configuration": {"endian": "little"},"name": "bytes"}, {"configuration": {"cname": "' + cname + '", "clevel": ' + str(
            clevel) + ', "shuffle": "' + shuffle + '", "typesize": ' + str(
            typesize) + ', "blocksize": ' + str(blocksize) + '},"name": "blosc"}]'
        print(codecs)
        self.pi2_write(codecs=codecs, chunk_shape=chunk_shape, data=data)
        read_arr = self.zarrita_read("test.zarr")
        assert np.array_equal(data, read_arr), "read_arr:\n " + str(read_arr) + " \n\narr:\n " + str(data)

    @pytest.mark.parametrize("chunk_shape", [[5, 5, 2]])
    @pytest.mark.parametrize("cname", ["lz4"])
    @pytest.mark.parametrize("clevel", [4])
    @pytest.mark.parametrize("shuffle", ["shuffle"])
    @pytest.mark.parametrize("typesize", [4])
    @pytest.mark.parametrize("blocksize", [0])
    @pytest.mark.parametrize("data", [full_arr, np.zeros_like(full_arr)])
    # @pytest.mark.blosc
    # TODO @pytest.mark.parametrize("cname", ["lz4", "lz4hc", "blosclz", "zstd", "snappy", "zlib"])
    def test_read_write_blosc(self, chunk_shape, cname, clevel, shuffle, typesize, blocksize, data):
        codecs = '[{"configuration": {"endian": "little"},"name": "bytes"}, {"configuration": {"cname": "' + cname + '", "clevel": ' + str(
            clevel) + ', "shuffle": "' + shuffle + '", "blocksize": ' + str(blocksize) + ', "typesize": ' + str(
            typesize) + '},"name": "blosc"}]'
        print(codecs)
        self.pi2_write(codecs=codecs, chunk_shape=chunk_shape, data=data)
        read_arr = self.pi2_read("test.zarr")
        assert np.array_equal(data, read_arr), "read_arr:\n " + str(read_arr) + " \n\narr:\n " + str(data)

    @pytest.mark.parametrize("separator", [".", "/", "-"])
    def test_read_separator(self, separator):
        self.zarrita_write(separator=separator)
        read_arr = self.pi2_read("zarrita.zarr")
        assert np.array_equal(self.arr, read_arr), "read_arr:\n " + str(read_arr) + " \n\narr:\n " + str(self.arr)

    blosc_codec = ',{"configuration": {"cname": "lz4", "clevel": 4, "shuffle": "shuffle", "blocksize": 0, "typesize": 4},"name": "blosc"}'

    @pytest.mark.parametrize("inner_chunk_shape", [[w, h, 1], [w, 1, 1], [1, 1, 1]])
    @pytest.mark.parametrize("shard_shape", [[w, h, d], [w, h, 1]])
    @pytest.mark.parametrize("bytes_bytes_codecs", ["", blosc_codec])
    @pytest.mark.parametrize("index_location", ["start", "end"])
    @pytest.mark.parametrize("data", [arr, arr0, arr42])
    def test_write_sharding(self, inner_chunk_shape, shard_shape, bytes_bytes_codecs, index_location, data, request):
        filename = request.node.name.replace('"', '') + ".zarr"
        config = '{"chunk_shape":' + str(
            inner_chunk_shape) + ',"codecs":[{"configuration":{"endian":"little"},"name":"bytes"}' + bytes_bytes_codecs + '],"index_codecs":[{"configuration":{"endian":"little"},"name":"bytes"}],"index_location":"' + index_location + '"}'
        codecs = '[{"name": "sharding_indexed", "configuration":' + config + '}]'
        print(codecs)
        self.pi2_write(codecs=codecs, chunk_shape=shard_shape, data=data, name=filename)
        read_arr = self.zarrita_read(filename)
        assert np.array_equal(data, read_arr), "read_arr:\n " + str(read_arr) + " \n\ndata:\n " + str(data)

    @pytest.mark.parametrize("inner_chunk_shape", [[w, h, 1], [w, 1, 1], [1, 1, 1]])
    @pytest.mark.parametrize("shard_shape", [[w, h, d], [w, h, 1]])
    @pytest.mark.parametrize("bytes_bytes_codecs", [[], [zarrita.codecs.blosc_codec(typesize=4)]])
    @pytest.mark.parametrize("index_location", ["start", "end"])
    @pytest.mark.parametrize("data", [arr, arr0, arr42])
    def test_read_sharding(self, inner_chunk_shape, shard_shape, bytes_bytes_codecs, index_location, data, request):
        filename = request.node.name.replace('"', '') + ".zarr"
        print(filename)
        self.zarrita_write(codecs=[
            zarrita.codecs.sharding_codec(
                chunk_shape=inner_chunk_shape,
                codecs=[
                    zarrita.codecs.bytes_codec(),
                    *bytes_bytes_codecs
                ],
                index_codecs=[zarrita.codecs.bytes_codec()],
                index_location=index_location
            ),
        ], chunk_shape=shard_shape, data=data, name=filename)
        read_arr = self.pi2_read(filename)
        assert np.array_equal(data, read_arr), "read_arr:\n " + str(read_arr) + " \n\ndata:\n " + str(data)

    blosc_codec = ',{"configuration": {"cname": "lz4", "clevel": 4, "shuffle": "shuffle", "blocksize": 0, "typesize": 4},"name": "blosc"}'

    @pytest.mark.parametrize("inner_chunk_shape", [[w, h, 1], [w, 1, 1], [1, 1, 1]])
    @pytest.mark.parametrize("shard_shape", [[w, h, d], [w, h, 1]])
    @pytest.mark.parametrize("bytes_bytes_codecs", ["", blosc_codec])
    @pytest.mark.parametrize("index_location", ["start", "end"])
    @pytest.mark.parametrize("data", [arr, arr0, arr42])
    def test_read_write_sharding(self, inner_chunk_shape, shard_shape, bytes_bytes_codecs, index_location, data, request):
        filename = request.node.name.replace('"', '') + ".zarr"
        config = '{"chunk_shape":' + str(
            inner_chunk_shape) + ',"codecs":[{"configuration":{"endian":"little"},"name":"bytes"}' + bytes_bytes_codecs + '],"index_codecs":[{"configuration":{"endian":"little"},"name":"bytes"}],"index_location":"' + index_location + '"}'
        codecs = '[{"name": "sharding_indexed", "configuration":' + config + '}]'
        print(codecs)
        self.pi2_write(codecs=codecs, chunk_shape=shard_shape, data=data, name=filename)
        read_arr = self.pi2_read(filename)
        assert np.array_equal(data, read_arr), "read_arr:\n " + str(read_arr) + " \n\ndata:\n " + str(data)

    def test_fill_value_in_api(self):
        pi2 = self.pi2
        name = "test_fill_value_in_api"
        data = np.zeros((2, 2, 2), dtype=np.float32)
        fill_value = 4
        shutil.rmtree(output_file(name), ignore_errors=True)
        write_img = pi2.newimage(pi2py2.ImageDataType.UInt32, list(data.shape))
        write_img.set_data(data.transpose(1, 0, 2))
        pi2.writezarr(write_img, output_file(name), list(data.shape),
                      '[{"configuration": {"endian": "little"},"name": "bytes"}]', fill_value, "/")
        store = zarrita.LocalStore(output_file(name))
        a = zarrita.Array.open(store)
        a = a.resize((3, 3, 3))
        assert a[:][2, 2, 2] == fill_value

    def test_dummy(self):
        pass


@pytest.mark.skip()
class Test_basic:

    pi2 = Pi2(library_path='../x64/Release')
    

    def check_result(self, is_ok, msg):
        """
        Checks test result
        """

        assert is_ok, msg



    def calc_difference(self, img1, img2):
        """
        Determines maximum absolute difference between img1 and img2.
        """

        pi2 = self.pi2

        imgf = pi2.newimage()
        img2f = pi2.newimage()

        pi2.convert(img1, imgf, ImageDataType.FLOAT32)
        pi2.convert(img2, img2f, ImageDataType.FLOAT32)

        pi2.subtract(imgf, img2f)
        pi2.abs(imgf)
        M = pi2.maxval(imgf)

        return M



    def check_distribution_test_result(self, file1, file2, operation_name, compname, tolerance=0.00001, data_type=ImageDataType.UNKNOWN):
        """
        Helper function used to check results of distributed processing tests.
        """

        pi2 = self.pi2

        # Load both results and calculate maximum absolute difference
        img1 = pi2.read(file1, data_type)
        img2 = pi2.read(file2, data_type)

        self.check_result(np.array_equal(img1.dimensions(), img2.dimensions()), f"Found difference in results of {compname} processing while testing {operation_name}. Image dimensions are not the same {file1} and {file2}.")

        M = self.calc_difference(img1, img2)

        # Check that the difference is zero
        self.check_result(M <= tolerance, f"Found difference in results of {compname} processing while testing {operation_name}. Max absolute difference = {M}.")
    


    def check_difference_delaying(self, testname, script, resultname='result', tolerance=0.00001, maxmem=30, chunk_size=[30, 32, 33]):
        """
        Calculates given script using distributed processing with delaying enabled and disabled.
        Compares the results and reports errors.
        """

        pi2 = self.pi2

        outfile_normal = output_file(f"{testname}_distributed_normal")
        outfile_delayed = output_file(f"{testname}_distributed_delayed")

        pi2.distribute(Distributor.LOCAL)
        pi2.maxmemory(maxmem)
        pi2.chunksize(chunk_size)

        pi2.delaying(False)
        pi2.run_script(script)
        pi2.writeraw(resultname, outfile_normal)

        pi2.delaying(True)
        pi2.run_script(script)
        pi2.writeraw(resultname, outfile_delayed)

        pi2.distribute(Distributor.NONE)

        self.check_distribution_test_result(outfile_normal, outfile_delayed, testname, 'normal distributed and delayed distributed', tolerance)




    def check_difference_normal_distributed(self, opname, args, resultname='result', infile=input_file(), tolerance=0.00001, convert_to_type=ImageDataType.UNKNOWN, maxmem=15, chunk_size=[64, 64, 64], max_jobs=0, out_prefix="out"):
        """
        Calculates operation normally and using distributed processing.
        Calculates difference between the results of the two versions and prints a message if the results do not match.
        opname is name of the function to call, e.g. 'gaussfilter'.
        args is array of parameters to the function. The source image is loaded to variables 'img' and 'result',
        and name of variable containing the result is stored in 'resultname' argument.
        If 'resultname' is empty, variable 'result' is assumed to store the result.
        """

        pi2 = self.pi2

        print(f"Testing {opname}...")
        print('-' * (len(opname) + 8 + 3))

        argstr = '_'.join(str(e) for e in args)
        argstr = argstr.replace('*', '_') # Remove non-filename characters
        argstr = argstr.replace('\n', '')
        outfile_normal = output_file(f"{out_prefix}_{opname}_{argstr}_normal")
        outfile_distributed = output_file(f"{out_prefix}_{opname}_{argstr}_distributed")

    


        # Run the operation using distributed processing
        pi2.distribute(Distributor.LOCAL)
        pi2.maxmemory(maxmem)
        pi2.chunksize(chunk_size)
        pi2.maxjobs(max_jobs)

        img = pi2.read(infile)
        if convert_to_type != ImageDataType.UNKNOWN:
            pi2.convert(img, convert_to_type)

        result = pi2.newlike(img)
        pi2.set(result, img)

        temp1 = pi2.newimage(ImageDataType.FLOAT32)
        temp2 = pi2.newimage(ImageDataType.FLOAT32)

        args_tmp = args.copy()
        for n, i in enumerate(args_tmp):
          if i == 'img':
             args_tmp[n] = img
          elif i == resultname:
              args_tmp[n] = result
          elif i == 'temp1':
              args_tmp[n] = temp1
          elif i == 'temp2':
              args_tmp[n] = temp2

        pi2.run_command(opname, args_tmp)

        for n, i in enumerate(args):
          if i == resultname:
              pi2.writeraw(args_tmp[n], outfile_distributed)


        pi2.distribute(Distributor.NONE)



        # Run the operation without distributed processing
        img = pi2.read(infile)
        if convert_to_type != ImageDataType.UNKNOWN:
            pi2.convert(img, convert_to_type)

        result = pi2.newlike(img)
        pi2.set(result, img)

        temp1 = pi2.newimage(ImageDataType.FLOAT32)
        temp2 = pi2.newimage(ImageDataType.FLOAT32)

        args_tmp = args.copy()
        for n, i in enumerate(args_tmp):
          if i == 'img':
             args_tmp[n] = img
          elif i == resultname:
              args_tmp[n] = result
          elif i == 'temp1':
              args_tmp[n] = temp1
          elif i == 'temp2':
              args_tmp[n] = temp2

        pi2.run_command(opname, args_tmp)

        for n, i in enumerate(args):
          if i == resultname:
              pi2.writeraw(args_tmp[n], outfile_normal)



        # Compare results
        self.check_distribution_test_result(outfile_normal, outfile_distributed, opname, 'distributed and local', tolerance, convert_to_type)



    
    def create_particle_labels_test_data(self):
        """
        Creates initial data for particle label growing test.
        """

        pi2 = self.pi2

        img = pi2.read(input_file('complicated_particles_1_38x36x21.raw'))

        pi2.set(img, [7, 7, 10], 2)
        pi2.set(img, [20, 18, 10], 3)
        pi2.set(img, [25, 25, 10], 4)

        pi2.writeraw(img, output_file('complicated_particles_point_labels'))


    
    def trace_skeleton_test(self, max_jobs=0):
        """
        Tests skeleton tracing in normal and distributed mode.
        """

        pi2 = self.pi2

        # Normal mode
        skele = pi2.read(input_file("real_skele_200x200x200.raw"))
        vertices = pi2.newimage()
        edges = pi2.newimage()
        measurements = pi2.newimage()
        points = pi2.newimage()
        pi2.tracelineskeleton(skele, vertices, edges, measurements, points)

        pi2.writeraw(vertices, output_file('vertices'))
        pi2.writeraw(edges, output_file('edges'))
        pi2.writeraw(measurements, output_file('measurements'))
        pi2.writeraw(points, output_file('points'))

        skele = pi2.read(input_file("real_skele_200x200x200.raw"))
        pi2.fillskeleton(skele, vertices, edges, measurements, points, 1)
        pi2.writeraw(skele, output_file('filled_skeleton'))
        pi2.clear()




        # Distributed mode
        pi2.distribute(Distributor.LOCAL)
        pi2.maxmemory(5)
        pi2.chunksize([100, 100, 100])
        pi2.maxjobs(max_jobs)

        skele = pi2.read(input_file("real_skele_200x200x200.raw"))
        vertices = pi2.newimage()
        edges = pi2.newimage()
        measurements = pi2.newimage()
        points = pi2.newimage()
        pi2.tracelineskeleton(skele, vertices, edges, measurements, points)

        pi2.writeraw(vertices, output_file('vertices_distributed'))
        pi2.writeraw(edges, output_file('edges_distributed'))
        pi2.writeraw(measurements, output_file('measurements_distributed'))
        pi2.writeraw(points, output_file('points_distributed'))

        skele = pi2.read(input_file("real_skele_200x200x200.raw"))
        pi2.fillskeleton(skele, vertices, edges, measurements, points, 1)
        pi2.writeraw(skele, output_file('filled_skeleton_distributed'))
        pi2.clear()

        pi2.distribute(Distributor.NONE)

        self.check_distribution_test_result(output_file('vertices'), output_file('vertices_distributed'), "tracelineskeleton", "vertices")
        # NOTE: Only vertices can be compared; the edge order might be different and probably the edges might be a bit different, too,
        # although bot are valid representations of the skeleton.
        #check_distribution_test_result(output_file('edges'), output_file('edges_distributed'), "tracelineskeleton", "edges")
        #check_distribution_test_result(output_file('measurements'), output_file('measurements_distributed'), "tracelineskeleton", "measurements")
        #check_distribution_test_result(output_file('points'), output_file('points_distributed'), "tracelineskeleton", "points", data_type='int32')
        # Skeleton filling should give the same result
        self.check_distribution_test_result(output_file('filled_skeleton'), output_file('filled_skeleton_distributed'), "tracelineskeleton", "filled skeleton")


    def test_trace_skeleton_multi_job(self):
        self.trace_skeleton_test(max_jobs=0)

    def test_trace_skeleton_one_job(self):
        self.trace_skeleton_test(max_jobs=1)




    def test_fill_skeleton(self):
        """
        Tests distributed and normal skeleton filling.
        """
        
        pi2 = self.pi2

        outfile_normal = output_file('filled_skeleton_length_normal')
        outfile_distributed = output_file('filled_skeleton_length_distributed')

        # First trace the skeleton.
        skele = pi2.read(input_file("real_skele_200x200x200.raw"))
        vertices = pi2.newimage()
        edges = pi2.newimage()
        measurements = pi2.newimage()
        points = pi2.newimage()
        pi2.tracelineskeleton(skele, vertices, edges, measurements, points)
    
        pi2.writeraw(vertices, output_file('vertices'))
        pi2.writeraw(edges, output_file('edges'))
        pi2.writeraw(measurements, output_file('measurements'))
        pi2.writeraw(points, output_file('points'))

        # Fill in normal mode
        skele = pi2.read(input_file("real_skele_200x200x200.raw"))
        pi2.fillskeleton(skele, vertices, edges, measurements, points, 1)
        pi2.writeraw(skele, outfile_normal)
        pi2.clear()

        # Fill in distributed mode
        pi2.distribute(Distributor.LOCAL)
        pi2.maxmemory(5)

        vertices = pi2.read(output_file('vertices'))
        edges = pi2.read(output_file('edges'))
        measurements = pi2.read(output_file('measurements'))
        points = pi2.read(output_file('points'), 'int32') # NOTE: We need to specify the data type here as otherwise it defaults to uint32. This is required only for .raw files.
        skele = pi2.read(input_file("real_skele_200x200x200.raw"))
        pi2.fillskeleton(skele, vertices, edges, measurements, points, 1)
        pi2.writeraw(skele, outfile_distributed)

        pi2.distribute(Distributor.NONE)

        self.check_distribution_test_result(outfile_normal, outfile_distributed, 'fillskeleton', 'distributed and local', 0)



    def test_morphorec(self):
        """
        Tests normal and distributed morphological reconstruction.
        """

        pi2 = self.pi2

        geom = pi2.read(input_file('complicated_particles_1_38x36x21.raw'))
    
        img = pi2.newlike(geom)
        pi2.set(img, [7, 7, 10], 2)
        pi2.set(img, [20, 18, 10], 3)
        pi2.set(img, [25, 25, 10], 4)
        pi2.writeraw(img, output_file('morphorec_seeds'))

        pi2.morphorec(img, geom)
        pi2.writeraw(img, output_file('morphorec_normal'))

        pi2.clear()
        pi2.distribute(Distributor.LOCAL)
        pi2.maxmemory(0.025)
        geom = pi2.read(input_file('complicated_particles_1_38x36x21.raw'))
        img = pi2.read(output_file('morphorec_seeds'))
        pi2.morphorec(img, geom)
        pi2.writeraw(img, output_file('morphorec_distributed'))

        pi2.distribute(Distributor.NONE)

        self.check_distribution_test_result(output_file('morphorec_normal'), output_file('morphorec_distributed'), 'morphorec', 'distributed and local', 0)



    
    def multimax_test(self, direction):
        """
        Tests normal and distributed max projection using second image.
        """

        pi2 = self.pi2

        img = pi2.read(input_file())
        imgVal = pi2.read(input_file_bin())

        maxproj = pi2.newimage(img.get_data_type())
        maxval = pi2.newimage(imgVal.get_data_type())
        pi2.maxproject(img, maxproj, direction, imgVal, maxval)

        pi2.writeraw(maxproj, output_file(f"max_projection_{direction}_normal"))
        pi2.writeraw(maxval, output_file(f"max_projection_value_{direction}_normal"))

        pi2.distribute(Distributor.LOCAL)
        pi2.maxmemory(15)

        img = pi2.read(input_file())
        imgVal = pi2.read(input_file_bin())

        maxproj = pi2.newimage(img.get_data_type())
        maxval = pi2.newimage(imgVal.get_data_type())
        pi2.maxproject(img, maxproj, direction, imgVal, maxval)

        pi2.writeraw(maxproj, output_file(f"max_projection_{direction}_distributed"))
        pi2.writeraw(maxval, output_file(f"max_projection_value_{direction}_distributed"))

        pi2.distribute(Distributor.NONE)

        self.check_distribution_test_result(output_file(f"max_projection_{direction}_normal"), output_file(f"max_projection_{direction}_distributed"), "multimax projection", "distributed and local", 0)
        self.check_distribution_test_result(output_file(f"max_projection_value_{direction}_normal"), output_file(f"max_projection_value_{direction}_distributed"), "multimax value", "distributed and local", 0)

    def test_multimax_x(self):
        self.multimax_test(0)

    def test_multimax_y(self):
        self.multimax_test(1)

    def test_multimax_z(self):
        self.multimax_test(2)






    def generate_particles(self):
        """
        analyze_particles function uses this to generate an image containing some random particles.
        """

        pi2 = self.pi2

        if pi2.isimagefile(output_file('particles')):
            return

        import random

        pi2.echo(False, False)

        img = pi2.newimage(ImageDataType.UINT16, 500, 500, 500)

        # Generate some particles

        # Spheres
        #for i in range(0, 1000):
        #    pos = [random.randint(0, 500), random.randint(0, 500), random.randint(0, 500)]
        #    r = random.randint(1, 20)
        #    col = random.randint(1, 255)
        #    pi2.sphere(img, pos, r, col)

        # Generic ellipsoids
        for i in range(0, 1000):
            pos = [random.randint(0, 500), random.randint(0, 500), random.randint(0, 500)]
            r = [random.randint(1, 20), random.randint(1, 20), random.randint(1, 20)]
            col = random.randint(1, 255)

            u1 = [random.randint(0, 500), random.randint(0, 500), random.randint(0, 500)]
            tmp = [random.randint(0, 500), random.randint(0, 500), random.randint(0, 500)]
            u2 = np.cross(u1, tmp)

            pi2.ellipsoid(img, pos, r, col, u1, u2)

        # Axis-aligned boxes
        #for i in range(0, 1000):
        #    pos = [random.randint(0, 500), random.randint(0, 500), random.randint(0, 500)]
        #    size = [random.randint(1, 20), random.randint(1, 20), random.randint(1, 20)]
        #    col = random.randint(1, 255)
        #    pi2.box(img, pos, size, col)

        # Generic boxes
        for i in range(0, 1000):
            pos = [random.randint(0, 500), random.randint(0, 500), random.randint(0, 500)]
            r = [random.randint(1, 20), random.randint(1, 20), random.randint(1, 20)]
            col = random.randint(1, 255)

            u1 = [random.randint(0, 500), random.randint(0, 500), random.randint(0, 500)]
            tmp = [random.randint(0, 500), random.randint(0, 500), random.randint(0, 500)]
            u2 = np.cross(u1, tmp)

            pi2.box(img, pos, r, col, u1, u2)


        pi2.writeraw(img, output_file('particles'))
        


    def analyze_particles_local(self, analyzers):
        """
        Helper for analyze_particles demo.
        """

        pi2 = self.pi2

        # Read generated data file
        img = pi2.read(output_file('particles'))
        pi2.threshold(img, 0)

        # Analyze particles
        result = pi2.newimage(ImageDataType.FLOAT32)
        pi2.analyzeparticles(img, result, analyzers)

        # Draw ellipsoids
        img = pi2.read(output_file('particles'))
        pi2.threshold(img, 0)
        pi2.drawellipsoids(img, analyzers, result, 2)
        pi2.writeraw(img, output_file('ellipsoids_local'))

        # Get a column to check that get_column works.
        test_column = pi2.get_column('Volume', result, analyzers)

        pa = result.get_data()
    

        return pa

    def analyze_particles_distributed(self, analyzers, max_jobs):
        """
        Helper for analyze_particles demo.
        """

        pi2 = self.pi2

        pi2.distribute(Distributor.LOCAL)
        pi2.maxmemory(50)
        pi2.maxjobs(max_jobs)

        # Read generated data file
        img = pi2.read(output_file('particles'))
        pi2.threshold(img, 0)

        # Analyze particles
        result = pi2.newimage(ImageDataType.FLOAT32)
        pi2.analyzeparticles(img, result, analyzers)

         # Draw ellipsoids
        img = pi2.read(output_file('particles'))
        pi2.threshold(img, 0)
        pi2.drawellipsoids(img, analyzers, result, 2)
        pi2.writeraw(img, output_file('ellipsoids_distr'))

        pa = result.get_data()

        pi2.distribute(Distributor.NONE)

        return pa


    def check_column(self, a, b, n, msg, tol):

        max_diff = 0
        for i in range(0, a.shape[0]):
            if np.isnan(a[i, n]) and np.isnan(b[i, n]):
                # Both nan, OK
                pass
            elif np.isposinf(a[i, n]) and np.isposinf(b[i, n]):
                # Both +inf, OK
                pass
            elif np.isneginf(a[i, n]) and np.isneginf(b[i, n]):
                # Both -inf, OK
                pass
            elif np.isnan(a[i, n]) and abs(b[i, n]) < 0.001:
                # OK: Allow one nan and one zero as those occur due to numerics.
                pass
            elif abs(a[i, n]) < 0.001 and np.isnan(b[i, n]):
                # OK: Allow one nan and one zero as those occur due to numerics.
                pass
            elif np.isnan(a[i, n]) and not np.isnan(b[i, n]):
                # One nan and one not, fail
                max_diff = np.nan
            elif not np.isnan(a[i, n]) and np.isnan(b[i, n]):
                # One nan and one not, fail
                max_diff = np.nan
            else:
                diff = np.abs(a[i, n] - b[i, n])
                if diff > max_diff:
                    max_diff = diff


        self.check_result((not np.isnan(diff)) and (max_diff <= tol), f"Local and distributed mode particle analysis differs in column {msg}, max absolute difference = {max_diff}")


    def analyze_particles(self, max_jobs=0):
        """
        Demonstrates particle (i.e, region or blob) analysis, and some drawing commands.
        """

        pi2 = self.pi2

        self.generate_particles()

        # Show analyzer names
        pi2.listanalyzers()

        # Analyze particles locally and using distributed computing mode
        analyzers = 'volume coordinates bounds boundingsphere isonedge pca'

        # Show titles of data columns
        print('Titles of columns in data table:')
        pi2.headers(analyzers)

        pa_local = self.analyze_particles_local(analyzers)

        # Flood fill all particles so that only background should be left
        img = pi2.read(output_file('particles'))
        pi2.threshold(img, 0)
        V = pa_local[:, 0]
        x = pa_local[:, 1]
        y = pa_local[:, 2]
        z = pa_local[:, 3]

        for i in range(0, len(x)):
            pi2.floodfill(img, [x[i], y[i], z[i]], 0)

        M = img.get_data().max()

        pi2.writeraw(img, output_file('black'))

        self.check_result(M == 0, "Not all particles are in the results table.")



        print(f"Local processing returned results array of size {pa_local.shape}")

        pa_dist = self.analyze_particles_distributed(analyzers, max_jobs)
        print(f"Distributed processing returned results array of size {pa_dist.shape}")


        # Particle order may be different in local and distributed processing, so sort both result arrays
        # before comparing.
        # NOTE: This will sort each column. That is not what we want!
        #pa_local.sort(0)
        #pa_dist.sort(0)
        # Sort according to volume and bounding box minx, miny, minz coordinates.
        pa_local = pa_local[np.lexsort((pa_local[:, 8], pa_local[:, 4], pa_local[:, 6], pa_local[:, 0]))]
        pa_dist = pa_dist[np.lexsort((pa_dist[:, 8], pa_dist[:, 4], pa_dist[:, 6], pa_dist[:, 0]))]

        np.savetxt(output_file('particle_analysis_local.csv'), pa_local)
        np.savetxt(output_file('particle_analysis_distributed.csv'), pa_dist)


        # Compare analysis results
    
        if not np.isclose(pa_local.shape, pa_dist.shape).all():
            self.check_result(False, 'Different number of results in local and distributed particle analysis.')
        else:
            # We need different tolerances for different columns!
            self.check_column(pa_local, pa_dist, 0, 'volume', 0)
            # Don't check X, Y, Z as those might differ as different points of the particle may be encountered first
            #self.check_column(pa_local, pa_dist, 1, 'X')
            #self.check_column(pa_local, pa_dist, 2, 'Y')
            #self.check_column(pa_local, pa_dist, 3, 'Z')
            self.check_column(pa_local, pa_dist, 4, 'minx', 0)
            self.check_column(pa_local, pa_dist, 5, 'maxx', 0)
            self.check_column(pa_local, pa_dist, 6, 'miny', 0)
            self.check_column(pa_local, pa_dist, 7, 'maxy', 0)
            self.check_column(pa_local, pa_dist, 8, 'minz', 0)
            self.check_column(pa_local, pa_dist, 9, 'maxz', 0)
            self.check_column(pa_local, pa_dist, 10, 'bounding sphere X', 0.01)
            self.check_column(pa_local, pa_dist, 11, 'bounding sphere Y', 0.01)
            self.check_column(pa_local, pa_dist, 12, 'bounding sphere Z', 0.01)
            self.check_column(pa_local, pa_dist, 13, 'bounding sphere radius', 0.01)
            self.check_column(pa_local, pa_dist, 14, 'is on edge', 0)
            self.check_column(pa_local, pa_dist, 15, 'CX', 0.01)
            self.check_column(pa_local, pa_dist, 16, 'CY', 0.01)
            self.check_column(pa_local, pa_dist, 17, 'CZ', 0.01)
            self.check_column(pa_local, pa_dist, 18, 'e', 0.01)
            self.check_column(pa_local, pa_dist, 19, 'l1', 0.01)
            self.check_column(pa_local, pa_dist, 20, 'l2', 0.01)
            self.check_column(pa_local, pa_dist, 21, 'l3', 0.01)
            # These angles are suprisingly inaccurate although the ellipsoids
            # are very nearly same (see the filling test)
            # 0.2 rad = 11.5 deg
            self.check_column(pa_local, pa_dist, 22, 'phi1', 0.2)
            self.check_column(pa_local, pa_dist, 23, 'theta1', 0.2)
            self.check_column(pa_local, pa_dist, 24, 'phi2', 0.2)
            self.check_column(pa_local, pa_dist, 25, 'theta2', 0.2)
            self.check_column(pa_local, pa_dist, 26, 'phi3', 0.2)
            self.check_column(pa_local, pa_dist, 27, 'theta3', 0.2)
            self.check_column(pa_local, pa_dist, 28, 'rmax', 0.01)
            self.check_column(pa_local, pa_dist, 29, 'd1', 1)
            self.check_column(pa_local, pa_dist, 30, 'd2', 1)
            self.check_column(pa_local, pa_dist, 31, 'd3', 1)
            self.check_column(pa_local, pa_dist, 32, 'bounding scale', 0.05)

        # Numerical precision causes small differences in the shapes of the ellipsoids, allow for that
        # This kind of expressions might be easier to write in NumPy as follows,
        # but then distributed processing possibility is of course lost:
        vis_local = pi2.read(output_file('ellipsoids_local')).get_data()
        vis_distr = pi2.read(output_file('ellipsoids_distr')).get_data()
        diff = np.abs(vis_local - vis_distr)
        diff_count = np.sum(diff > 0)
        self.check_result(diff_count < 100, "More than 100 pixels are different between ellipsoid visualizations made using local and distributed processing.")



    def test_analyze_particles_multi_job(self):
        self.analyze_particles(max_jobs=0)

    def test_analyze_particles_single_job(self):
        self.analyze_particles(max_jobs=1)



    def test_analyze_labels(self):

        pi2 = self.pi2

        self.generate_particles()

        img = pi2.read(output_file('particles'))
        analyzers = 'volume coordinates'
        pa_local_img = pi2.newimage(ImageDataType.FLOAT32)
        pi2.analyzelabels(img, pa_local_img, analyzers)
        pa_local = pa_local_img.get_data()

        # Fill all particles so that only background should be left
        V = pa_local[:, 0]
        x = pa_local[:, 1]
        y = pa_local[:, 2]
        z = pa_local[:, 3]

        data = img.get_data()

        for i in range(0, len(x)):
            pix = data[int(y[i]), int(x[i]), int(z[i])] # NOTE: the order of coordinates is y, x, z in the Numpy array!
            pi2.replace(img, pix, 0)

        M = img.get_data().max()

        pi2.writeraw(img, output_file('black'))

        self.check_result(M == 0, "Not all particles are in the results table.")


    def test_dimension_broadcast(self):

        pi2 = self.pi2

        img = pi2.read(input_file())
        img2 = pi2.newlike(img, width=img.get_width(), height=img.get_height(), depth=1)
        pi2.set(img2, 2)
        pi2.divide(img, img2, True)

        pi2.writeraw(img, output_file('broadcasted_divide'))

        img_gt = pi2.read(input_file())
        pi2.divide(img_gt, 2)
        pi2.writeraw(img, output_file('normal_divide'))
    
        pi2.distribute(Distributor.LOCAL)
        pi2.maxmemory(1)

        img = pi2.read(input_file())
        img2 = pi2.newlike(img, width=img.get_width(), height=img.get_height(), depth=1)
        pi2.set(img2, 2)
        pi2.divide(img, img2, True)

        pi2.writeraw(img, output_file('broadcasted_distributed_divide'))

        pi2.distribute(Distributor.NONE)

        self.check_distribution_test_result(output_file('broadcasted_divide'), output_file('normal_divide'), "broadcasted division", "broadcasted to normal")
        self.check_distribution_test_result(output_file('broadcasted_distributed_divide'), output_file('normal_divide'), "broadcasted division", "distributed to normal")


    def test_rotate(self):

        pi2 = self.pi2

        img = pi2.read(input_file())
        out = pi2.newimage()
        pi2.rot90cw(img, out)
        pi2.writeraw(out, output_file('rot90cw'))

        # TODO: Asserts?



    def test_twoimage_distribution(self):

        pi2 = self.pi2

        img = pi2.read(input_file())
        out = pi2.newimage()

        pi2.finalizetmap(img, out)

        pi2.writeraw(out, output_file('finalizetmap_normal'))

        pi2.distribute(Distributor.LOCAL)
        pi2.maxmemory(15)

        img = pi2.read(input_file())
        out = pi2.newimage()

        pi2.finalizetmap(img, out)

        pi2.writeraw(out, output_file('finalizetmap_distributed'))
        pi2.distribute(Distributor.NONE)

        self.check_distribution_test_result(output_file('finalizetmap_normal'), output_file('finalizetmap_distributed'), 'two-image input-output commands distribution', 'distributed and local')



    def test_histogram(self):

        pi2 = self.pi2

        img = pi2.newimage(ImageDataType.FLOAT32, 200, 200, 200)
        pi2.ramp(img, 0)
        pi2.noise(img, 0, 10000000)
        pi2.writeraw(img, output_file('hist_input'))


        pi2.distribute(Distributor.LOCAL)
        pi2.maxmemory(15)
    
        img = pi2.read(output_file('hist_input'))

        hst1 = pi2.newimage(ImageDataType.FLOAT32)
        bins1 = pi2.newimage(ImageDataType.FLOAT32)
        pi2.hist(img, hst1, bins1, 0.0, 1000.0, 500)

        hst2 = pi2.newimage(ImageDataType.FLOAT32)
        bins2 = pi2.newimage(ImageDataType.FLOAT32)
        pi2.hist(img, hst2, bins2, 0.0, 1000.0, 250)

        hdata1 = hst1.get_data()
        hdata2 = hst2.get_data()

        print(hdata1.sum())
        print(hdata2.sum())

        self.check_result(hdata1.sum() == hdata2.sum(), "Histogram normalization does not match.")
        
        pi2.distribute(Distributor.NONE)




    def test_autothreshold(self):

        pi2 = self.pi2

        def single_test(method):

            img = pi2.read(input_file())
            pi2.linmap(img, 0, 600, 0, 255)
            pi2.convert(img, ImageDataType.UINT8)
            pi2.writetif(img, output_file('threshold_orig'))

            th = pi2.newlike(img)
            pi2.set(th, img)
            pi2.autothreshold(th, method)
            pi2.multiply(th, 255)
            pi2.writetif(th, output_file(f"threshold_{method}"))

            # TODO: Asserts?

        single_test(AutoThresholdMethod.OTSU)
        single_test(AutoThresholdMethod.HUANG)
        single_test(AutoThresholdMethod.INTERMODES)
        single_test(AutoThresholdMethod.ISODATA)
        single_test(AutoThresholdMethod.LI)
        single_test(AutoThresholdMethod.MAXENTROPY)
        single_test(AutoThresholdMethod.MEAN)
        single_test(AutoThresholdMethod.MINERROR)
        single_test(AutoThresholdMethod.MINIMUM)
        single_test(AutoThresholdMethod.MOMENTS)
        single_test(AutoThresholdMethod.PERCENTILE)
        single_test(AutoThresholdMethod.RENYI)
        single_test(AutoThresholdMethod.SHANBHAG)
        single_test(AutoThresholdMethod.TRIANGLE)
        single_test(AutoThresholdMethod.YEN)
        single_test(AutoThresholdMethod.MEDIAN)
        single_test(AutoThresholdMethod.MIDGREY)
        single_test(AutoThresholdMethod.NIBLACK)
        single_test(AutoThresholdMethod.PHANSALKAR)
        single_test(AutoThresholdMethod.SAUVOLA)
        single_test(AutoThresholdMethod.BERNSEN)



    def test_tif_and_tiff(self):
        """
        Ensures that both .tif and .tiff sequences can be read.
        """

        pi2 = self.pi2

        img = pi2.read(input_file())

        pi2.writesequence(img, output_file("./sequence_tif/img_@.tif"))
        pi2.writesequence(img, output_file("./sequence_tiff/img_@.tiff"))

        img2 = pi2.read(output_file("./sequence_tif/img_@.tif"))
        img3 = pi2.read(output_file("./sequence_tiff/img_@.tiff"))

        M1 = self.calc_difference(img, img2)
        M2 = self.calc_difference(img, img3)

        # Check that the difference is zero
        self.check_result(M1 <= 0, f"ERROR: Difference in tif and tiff sequence reading.")
        self.check_result(M2 <= 0, f"ERROR: Difference in tif and tiff sequence reading.")


    def get_pixels(self, max_jobs=0):
        """
        Checks that various ways of getting pixels from the (normal and distributed) images give the same result.
        """

        pi2 = self.pi2

        pi2.distribute(Distributor.NONE)
        # Generate image
        w = 100
        h = 80
        d = 120
        img = pi2.newimage(ImageDataType.UINT16, w, h, d)
        pi2.ramp(img, 0)
        pi2.noise(img, 0, 25)
        pi2.writeraw(img, output_file("ramp"))

        # Generate some positions
        N = 300
        #positions = pi2.newimage(ImageDataType.FLOAT32, 3, N)
        #pos = positions.get_data()
        pos = np.zeros([3, N])
        for i in range(0, N):
            pos[0, i] = random.randint(0, w-2) + 0.5
            pos[1, i] = random.randint(0, h-1)
            pos[2, i] = random.randint(0, d-1)
        #pos = np.array([[34.5, 42, 13]])

        # Get pixels at positions (non-distributed)
        img = pi2.read(output_file("ramp"))
        out = pi2.newimage(img.get_data_type())
        pi2.get(img, out, pos)
        data_normal = out.get_data()


        # Get pixels at positions (distributed)
        pi2.distribute(Distributor.LOCAL)
        pi2.maxmemory(0.25)
        pi2.maxjobs(max_jobs)

        img = pi2.read(output_file("ramp"))
        out = pi2.newimage(img.get_data_type())
        pi2.get(img, out, pos)
        data_distributed = out.get_data()

        pi2.distribute(Distributor.NONE)

        # Get pixels at positions (through NumPy)
        img = pi2.read(output_file("ramp"))
        pyimg = img.to_numpy()

        N = pos.shape[1]
        data_numpy = np.zeros(N)
        for i in range(0, N):
            data_numpy[i] = pyimg[int(pos[0, i] + 0.5), int(pos[1, i] + 0.5), int(pos[2, i] + 0.5)]

        self.check_result(np.isclose(data_normal, data_distributed).all(), "get pixel normal != distributed")
        self.check_result(np.isclose(data_normal, data_numpy).all(), "get pixel normal != numpy")

        print("point: normal = numpy = distributed")
        for i in range(0, N):
            print(f"{pos[:, i]}: {data_normal[i]} = {data_numpy[i]} = {data_distributed[i]}")


    def test_get_pixels_multi_job(self):
        self.get_pixels(max_jobs=0)

    def test_get_pixels_single_job(self):
        self.get_pixels(max_jobs=1)


    def test_set_pixels(self):

        pi2 = self.pi2

        img = pi2.newimage(ImageDataType.UINT8, 100, 200, 300)
        p = np.array([[0, 0, 0], [10, 20, 30], [70, 70, 70]], dtype=float)
        v = np.array([[10, 20, 30]], dtype=np.uint8)

        pos = pi2.newimage()
        pos.set_data(p)

        gt = pi2.newimage()
        gt.set_data(v)

        pi2.set(img, pos, gt)

        pi2.writetif(img, output_file("set_pixels"))


        v2 = pi2.newimage()
        pi2.get(img, v2, pos)

    

        M = self.calc_difference(gt, v2)

        self.check_result(M <= 0, "Pixels are not written and read correctly.")




    def test_distributed_numpy(self):
        """
        Writing numpy arrays to disk caused exception in distributed mode.
        """

        pi2 = self.pi2

        pi2.distribute(Distributor.LOCAL)

        nparr = np.zeros([100, 100, 100])

        pi2.writeraw(nparr, output_file("np_distributed"))

        pi2.distribute(Distributor.NONE)

        # TODO: Asserts


    def test_memory(self):
        """
        Checks that image memory is freed when variables are cleared.
        """

        # NOTE: This is not 100 % reliable test for some reason. Sometimes it fails for no reason.

        pi2 = self.pi2

        import os
        import psutil
        process = psutil.Process(os.getpid())

        mem = process.memory_info().rss
        self.check_result(mem < 250 * 1024 * 1024, "consuming too much memory before allocation of an image")

        img = pi2.newimage(ImageDataType.UINT8, 1000, 1000, 1000)

        mem = process.memory_info().rss
        self.check_result(mem > 500 * 1024 * 1024, "consuming too little memory after allocation of an image")

        pi2.clear(img.name)

        mem = process.memory_info().rss
        self.check_result(mem < 250 * 1024 * 1024, "consuming too much memory after de-allocation of an image")



    def test_metadata(self):
        """
        Tests storing and retrieving image meta-data.
        """

        pi2 = self.pi2

        img = pi2.newimage(ImageDataType.UINT8, 1, 1, 1)
        str = pi2.newvalue("string")

        pi2.setmeta(img, "key1", "value1")
        pi2.setmeta(img, "key2", "value2")
        pi2.getmeta(img, "key1", str)

        self.check_result(str.as_string() == "value1", "value")
    
        pi2.listmeta(img, str)
        self.check_result(str.as_string() == "key1, key2", "listmeta before")

        pi2.writemeta(img, output_file("metadata"))
        pi2.clearmeta(img)

        pi2.listmeta(img, str)
        self.check_result(str.as_string() == "", "meta key list not empty")

        pi2.readmeta(img, output_file("metadata"))
        pi2.listmeta(img, str)
        self.check_result(str.as_string() == "key1, key2", "meta key list contains incorrect values")

        pi2.getmeta(img, "key1", str)
        self.check_result(str.as_string() == "value1", "key1 after read")

        pi2.getmeta(img, "key2", str)
        self.check_result(str.as_string() == "value2", "key2 after read")



    def test_named_variables(self):
        """
        Tests named non-image variables.
        """

        pi2 = self.pi2

        img = pi2.newimage(ImageDataType.UINT8, 1, 1, 1)
        str_key = pi2.newvalue("string")
        str_value = pi2.newvalue("string", "VALUE")

        self.check_result(str_key.as_string() == "", "key before doing anything")
        self.check_result(str_value.as_string() == "VALUE", "value before doing anything")

        pi2.set(str_key, "KEY")

        self.check_result(str_key.as_string() == "KEY", "key after setting")

        pi2.setmeta(img, str_key, str_value)

        self.check_result(str_key.as_string() == "KEY", "key after setmeta")
        self.check_result(str_value.as_string() == "VALUE", "value after setmeta")

        pi2.set(str_value, "---")
        self.check_result(str_value.as_string() == "---", "value after set")

        pi2.getmeta(img, str_key, str_value)
    
        self.check_result(str_key.as_string() == "KEY", "key after getmeta")
        self.check_result(str_value.as_string() == "VALUE", "value after getmeta")


    def test_set_overloads(self):

        pi2 = self.pi2

        img1 = pi2.newimage(ImageDataType.UINT16)
        img2 = pi2.newlike(img1)

        pi2.set(img2, img1)

        # TODO: Asserts?


    def test_big_tiff_write(self):
        pi2 = self.pi2

        img = pi2.newimage(ImageDataType.UINT8, 5*1024, 1024, 1024)
        pi2.ramp(img, 0)
        pi2.writetif(img, output_file("big_tiff.tif"))
        pi2.writeraw(img, output_file("big_raw"))

        # TODO: Asserts?

    def test_big_tiff(self):

        pi2 = self.pi2

        self.test_big_tiff_write()

        # Whoops: I don't have enough memory for this in my computer...
        img1 = pi2.read(output_file("big_tiff.tif"))
        img2 = pi2.read(output_file("big_raw"))
        M = self.calc_difference(img1, img2)
        self.check_result(M <= 0, f"ERROR: Images saved as .tif and .raw are not equal.")



    def test_lz4_files(self):

        pi2 = self.pi2

        img1 = pi2.newimage(ImageDataType.FLOAT32, 100, 250, 200)
        pi2.ramp3(img1)
        pi2.writelz4(img1, output_file("lz4test"))
        pi2.writetif(img1, output_file("lz4test.tif"))

        img2 = pi2.read(output_file("lz4test.lz4raw"))
        self.check_result(self.calc_difference(img1, img2) == 0, "Image saved and read from .lz4raw dataset changed in the I/O process.")
    

    def test_nn5_files(self):

        pi2 = self.pi2

        img1 = pi2.newimage(ImageDataType.FLOAT32, 100, 250, 200)
        pi2.ramp3(img1)
        pi2.writenn5(img1, output_file("nn5test"))
        pi2.writenn5(img1, output_file("nn5test_small_chunks"), [20, 40, 60])
        pi2.writetif(img1, output_file("nn5test.tif"))

        img2 = pi2.read(output_file("nn5test"))
        self.check_result(self.calc_difference(img1, img2) == 0, "Image saved and read from NN5 dataset changed in the I/O process.")

        img3 = pi2.read(output_file("nn5test_small_chunks"))
        self.check_result(self.calc_difference(img1, img3) == 0, "Image saved and read from NN5 small chunks dataset changed in the I/O process.")


    def test_dead_pixels(self):

        pi2 = self.pi2

        img = pi2.newimage(ImageDataType.FLOAT32, 100, 100)
        x = [20, 20, 0]
        pi2.set(img, x, float('nan'))
        pi2.deadpixelremoval(img)
        val = img.get_data()[x[0], x[1]]
        self.check_result(val == 0, "dead pixel not removed")



    def test_lsf_cluster(self):

        if os.name == 'nt':
            pytest.skip("Unable to test LSF cluster functionality in Windows.")

        pi2 = self.pi2

        pi2.distribute(Distributor.LSF)
        pi2.maxmemory(1)
        img = pi2.newimage(ImageDataType.UINT16, 100, 100, 100)
        pi2.add(img, 10)
        pi2.writeraw(img, output_file("lsf/result"))

        pi2.distribute(Distributor.NONE)


    def test_max_filters_dist(self):
        for r in range(1, 10):
            self.check_difference_normal_distributed('maxfilter', ['img', 'result', r, True, 'ellipsoidal', 'zero'])
            self.check_difference_normal_distributed('maxfilter', ['img', 'result', r, True, 'ellipsoidal', 'nearest'])
            
    def test_min_filters_dist(self):
        for r in range(1, 10):
            self.check_difference_normal_distributed('minfilter', ['img', 'result', r, True, 'ellipsoidal', 'zero'])
            self.check_difference_normal_distributed('minfilter', ['img', 'result', r, True, 'ellipsoidal', 'nearest'])
       
    def test_opening_filters_dist(self):
        for r in range(1, 10):
            self.check_difference_normal_distributed('openingfilter', ['result', r, True, 'ellipsoidal', 'zero'])
            self.check_difference_normal_distributed('openingfilter', ['result', r, True, 'ellipsoidal', 'nearest'])
   
    def test_closing_filters_dist(self):
        for r in range(1, 10):
            self.check_difference_normal_distributed('closingfilter', ['result', r, True, 'ellipsoidal', 'zero'])
            self.check_difference_normal_distributed('closingfilter', ['result', r, True, 'ellipsoidal', 'nearest'])

    def test_floodfill_dist(self):
        self.check_difference_normal_distributed('floodfill', ['img', [0, 0, 0], 100], 'img', input_file_bin(), maxmem=5)
    
    def generate_floodfill_geometry(self):
        pi = self.pi2

        # Generate geometry
        img = pi.newimage(ImageDataType.UInt8, 200, 200, 200)
        pi.box(img, [10, 0, 0], [30, 150, 200], 128)
        pi.box(img, [10, 130, 0], [150, 20, 200], 128)
        geom_file = output_file('floodfill_geometry1')
        pi.writeraw(img, geom_file)
        return geom_file

    def generate_floodfill_geometry_closed(self):
        pi = self.pi2

        # Generate geometry
        img = pi.newimage(ImageDataType.UInt8, 200, 200, 200)
        pi.box(img, [10, 0, 0], [30, 150, 200], 128)
        pi.box(img, [10, 130, 0], [500, 20, 200], 128)
        pi.box(img, [10, 0, 0], [500, 150, 10], 128)
        geom_file = output_file('floodfill_geometry2')
        pi.writeraw(img, geom_file)
        return geom_file


    def test_floodfill_2_dist(self):
        geom_file = self.generate_floodfill_geometry()
        self.check_difference_normal_distributed('floodfill', ['img', [0, 0, 0], 200], 'img', geom_file, maxmem=2, chunk_size=[100, 100, 100], out_prefix="chunk_100")
        
    def test_floodfill_3_dist(self):
        geom_file = self.generate_floodfill_geometry()  
        self.check_difference_normal_distributed('floodfill', ['img', [0, 0, 0], 200], 'img', geom_file, maxmem=2, chunk_size=[64, 64, 64], out_prefix="chunk_64")
        
    def test_floodfill_4_dist(self):
        geom_file = self.generate_floodfill_geometry()
        self.check_difference_normal_distributed('floodfill', ['img', [0, 0, 0], 200], 'img', geom_file, maxmem=2, chunk_size=[60, 100, 200], out_prefix="chunk_mix")
        
    def test_floodfill_5_dist(self):
        geom_file = self.generate_floodfill_geometry_closed()
        self.check_difference_normal_distributed('floodfill', ['img', [0, 0, 0], 200], 'img', geom_file, maxmem=2, chunk_size=[60, 100, 200], out_prefix="chunk_mix_closed")


    def test_floodfill_multiseed_dist_1(self):
        geom_file = self.generate_floodfill_geometry()
        self.check_difference_normal_distributed('floodfill', ['img', np.array([[0, 0, 0], [199, 199, 199]], dtype=np.float32).transpose(), 200], 'img', geom_file, maxmem=4, chunk_size=[100, 100, 100], out_prefix="multiseed")

    def test_floodfill_multiseed_dist_2(self):
        geom_file = self.generate_floodfill_geometry()
        self.check_difference_normal_distributed('floodfill', ['img', np.array([[0, 0, 0], [12, 1, 1]], dtype=np.float32).transpose(), 200], 'img', geom_file, maxmem=4,  chunk_size=[100, 100, 100], out_prefix="multiseed")


    def test_floodfill_multiseed_dist_3(self):
        geom_file = self.generate_floodfill_geometry_closed()
        self.check_difference_normal_distributed('floodfill', ['img', np.array([[0, 0, 0], [199, 199, 199]], dtype=np.float32).transpose(), 200], 'img', geom_file, maxmem=4, chunk_size=[100, 100, 100], out_prefix="multiseed_closed")

    def test_floodfill_multiseed_dist_4(self):
        geom_file = self.generate_floodfill_geometry_closed()
        self.check_difference_normal_distributed('floodfill', ['img', np.array([[0, 0, 0], [12, 1, 1]], dtype=np.float32).transpose(), 200], 'img', geom_file, maxmem=4,  chunk_size=[100, 100, 100], out_prefix="multiseed_closed")


    def test_gauss_dist(self):
        self.check_difference_normal_distributed('gaussfilter', ['img', 'result', 2])
    
    def test_bin_dist(self):
        self.check_difference_normal_distributed('bin', ['img', 'result', 9])

    def test_scale_dist(self):
        self.check_difference_normal_distributed('scale', ['img', 'result', 0], tolerance=1)

    def test_scale1_dist(self):
        self.check_difference_normal_distributed('scale', ['img', 'result', 0.2, False, Interpolation.NEAREST], tolerance=1)
    
    def test_scale2_dist(self):
        self.check_difference_normal_distributed('scale', ['img', 'result', 0.2, True, Interpolation.NEAREST], tolerance=1, maxmem=35)
    
    def test_scale3_dist(self):
        self.check_difference_normal_distributed('scale', ['img', 'result', 0.2, False, Interpolation.LINEAR], tolerance=1)
    
    def test_scale4_dist(self):
        self.check_difference_normal_distributed('scale', ['img', 'result', 0.2, True, Interpolation.LINEAR], tolerance=1, maxmem=35)
    
    def test_scale5_dist(self):
        self.check_difference_normal_distributed('scale', ['img', 'result', 0.2, False, Interpolation.CUBIC], tolerance=1)
    
    def test_scale5_dist(self):
        self.check_difference_normal_distributed('scale', ['img', 'result', 0.2, True, Interpolation.CUBIC], tolerance=1, maxmem=35)
    
    def test_scale6_dist(self):
        self.check_difference_normal_distributed('scale', ['img', 'result', 0.7, False, Interpolation.NEAREST], tolerance=1)
    
    def test_scale7_dist(self):
        self.check_difference_normal_distributed('scale', ['img', 'result', 0.7, True, Interpolation.NEAREST], tolerance=1)
    
    def test_scale8_dist(self):
        self.check_difference_normal_distributed('scale', ['img', 'result', 0.7, False, Interpolation.LINEAR], tolerance=1)
    
    ## Failing tests:
    ## The failure has probably something to do with the processing block size being slightly incorrect if averaging
    ## when downscaling is used, if the processing block size happens to fulfill some condition. The results are
    ## very near to the correct ones, however.
    ##
    ##def test_scale9_dist(self):
    ##    self.check_difference_normal_distributed('scale', ['img', 'result', 0.7, True, Interpolation.LINEAR], tolerance=1)
    ##
    ##def test_scale10_dist(self):
    ##    self.check_difference_normal_distributed('scale', ['img', 'result', 0.7, True, Interpolation.CUBIC], tolerance=1)
    
    def test_scale11_dist(self):
        self.check_difference_normal_distributed('scale', ['img', 'result', 0.7, False, Interpolation.CUBIC], tolerance=1)
    
    def test_scale12_dist(self):
        self.check_difference_normal_distributed('scale', ['img', 'result', 2.1, False, Interpolation.NEAREST], tolerance=1)
    
    def test_scale13_dist(self):
        self.check_difference_normal_distributed('scale', ['img', 'result', 2.1, True, Interpolation.NEAREST], tolerance=1)
    
    def test_scale14_dist(self):
        self.check_difference_normal_distributed('scale', ['img', 'result', 2.1, False, Interpolation.LINEAR], tolerance=1)
    
    def test_scale15_dist(self):
        self.check_difference_normal_distributed('scale', ['img', 'result', 2.1, True, Interpolation.LINEAR], tolerance=1)
    
    def test_scale16_dist(self):
        self.check_difference_normal_distributed('scale', ['img', 'result', 2.1, False, Interpolation.CUBIC], tolerance=1)
    
    def test_scale17_dist(self):
        self.check_difference_normal_distributed('scale', ['img', 'result', 2.1, True, Interpolation.CUBIC], tolerance=1)
    
    def test_scalelabels1_dist(self):
        self.check_difference_normal_distributed('scalelabels', ['img', 'result', 2])
    
    def test_scalelabels2_dist(self):
        self.check_difference_normal_distributed('scalelabels', ['img', 'result', 3], maxmem=50)
    
    def test_crop_dist(self):
        self.check_difference_normal_distributed('crop', ['img', 'result', '[10, 20, 30]', '[110, 120, 130]'])
    
    def test_sumproject_dist(self):
        self.check_difference_normal_distributed('sumproject', ['img', 'result', 2])
    
    def test_meanproject_dist(self):
        self.check_difference_normal_distributed('meanproject', ['img', 'result', 1])
    
    def test_maxproject_dist(self):
        self.check_difference_normal_distributed('maxproject', ['img', 'result', 0])
    
    def test_sum1_dist(self):
        self.check_difference_normal_distributed('sum', ['img', 'result', True], tolerance=0)
    
    def test_sum2_dist(self):
        self.check_difference_normal_distributed('sum', ['img', 'result', True], tolerance=0.000001e8, convert_to_type=ImageDataType.FLOAT32)
    
    def test_mean_dist(self):
        self.check_difference_normal_distributed('mean', ['img', 'result'], tolerance=0.001)
    
    def test_maxval_dist(self):
        self.check_difference_normal_distributed('maxval', ['img', 'result', True])
    
    def test_satofilter_dist(self):
        self.check_difference_normal_distributed('satofilter', ['img', 'result', 1], maxmem=100)
    
    def test_dualthreshold_dist(self):
        self.check_difference_normal_distributed('dualthreshold', ['img', 300, 500], 'img')
    
    def test_regionremoval_dist(self):
        self.check_difference_normal_distributed('regionremoval', ['img', 500, 'all'], 'img', input_file_bin())
    
    def test_dmap2_dist(self):
        self.check_difference_normal_distributed('dmap2', ['img', 'result'], 'result', input_file_bin(), convert_to_type=ImageDataType.UINT16)
    
    def test_dmap_dist(self):
        self.check_difference_normal_distributed('dmap', ['img', 'result'], 'result', input_file_bin(), convert_to_type=ImageDataType.UINT16)
    
    def test_dmap2_2_dist(self):
        self.check_difference_normal_distributed('dmap2', ['img', 'result'], 'result', input_file_bin(), convert_to_type=ImageDataType.UINT32)
    
    def test_dmap_2_dist(self):
        self.check_difference_normal_distributed('dmap', ['img', 'result'], 'result', input_file_bin(), convert_to_type=ImageDataType.UINT32)
    
    def test_dmap2_3_dist(self):
        self.check_difference_normal_distributed('dmap2', ['img', 'result'], 'result', input_file_bin(), convert_to_type=ImageDataType.FLOAT32)
    
    def test_dmap_3_dist(self):
        self.check_difference_normal_distributed('dmap', ['img', 'result'], 'result', input_file_bin(), convert_to_type=ImageDataType.FLOAT32)
    
    def test_dmap_4_dist(self):
        self.check_difference_normal_distributed('dmap', ['img', 'img'], 'img', input_file_bin(), convert_to_type=ImageDataType.FLOAT32)
    
    def test_derivative1_dist(self):
        self.check_difference_normal_distributed('derivative', ['img', 'result', 1,  0], 'result', maxmem=20)
    
    def test_derivative2_dist(self):
        self.check_difference_normal_distributed('derivative', ['img', 'result', 1,  1], 'result', maxmem=20)
    
    def test_derivative3_dist(self):
        self.check_difference_normal_distributed('derivative', ['img', 'result', 1,  2], 'result', maxmem=20)
    
    def test_canny1_dist(self):
        self.check_difference_normal_distributed('canny', ['img', 1,  10, 100], 'img', maxmem=20, max_jobs=0)
    
    def test_canny2_dist(self):
        self.check_difference_normal_distributed('canny', ['img', 1,  10, 100], 'img', maxmem=20, max_jobs=1)
    
    def test_convert1_dist(self):
        self.check_difference_normal_distributed('convert', ['img', 'out', ImageDataType.FLOAT32], 'out')
    
    def test_convert2_dist(self):
        self.check_difference_normal_distributed('convert', ['img', 'img', ImageDataType.FLOAT32], 'img')
    
    def test_convert3_dist(self):
        self.check_difference_normal_distributed('convert', ['img', ImageDataType.FLOAT32], 'img')
    
    def test_hist1_dist(self):
        self.check_difference_normal_distributed('hist', ['img', 'result', 'temp1', 0, 1000, 1000], 'result')
    
    def test_hist2_dist(self):
        self.check_difference_normal_distributed('hist', ['img', 'result', 'temp1', 0, 500, 500], 'result')
    
    def test_whist_dist(self):
        self.check_difference_normal_distributed('whist', ['img', 'img', 'result', 'temp1', 0, 500, 500], 'result', convert_to_type=ImageDataType.FLOAT32)
    
    def test_hist2_2_dist(self):
        self.check_difference_normal_distributed('hist2', ['img', 0, 1000, 1000, 'img', 0, 100, 100, 'result', 'temp1', 'temp2'], 'result')
    
    def test_whist2_2_dist(self):
        self.check_difference_normal_distributed('whist2', ['img', 0, 1000, 1000, 'img', 0, 100, 100, 'img', 'result', 'temp1', 'temp2'], 'result', convert_to_type=ImageDataType.FLOAT32)
    
    def test_danielsson2_dist(self):
        self.check_difference_normal_distributed('danielsson2', ['img', 'result'], 'result', input_file('t1-head_bin_dmap_256x256x129.raw'))
    
    def test_maskedmean_dist(self):
        self.check_difference_normal_distributed('maskedmean', ['img', 'result', 0], 'result', input_file_bin(), convert_to_type=ImageDataType.FLOAT32)
    
    def test_tmap1_dist(self):
        self.check_difference_normal_distributed('tmap', ['img', 'result', 0, False, False], 'result', input_file_bin(), convert_to_type=ImageDataType.UINT16)
    
    def test_tmap2_dist(self):
        self.check_difference_normal_distributed('tmap', ['img', 'result', 0, False, False, '[50, 50, 50]'], 'result', input_file_bin(), convert_to_type=ImageDataType.UINT16)
    
    def test_growlabels_dist(self):
        self.create_particle_labels_test_data()
        self.check_difference_normal_distributed('growlabels', ['img', 1, 0], 'img', output_file('complicated_particles_point_labels'), maxmem=0.03)
    
    def test_cylinderorientation_dist(self):
        ## We do this test in two parts (we still check only theta or phi and not both, but if one is ok, the other should be ok, too!)
        self.check_difference_normal_distributed('cylinderorientation', ['img', 'result', 'result', 1, 1], 'result', input_file(), convert_to_type=ImageDataType.FLOAT32, maxmem=100)
        self.check_difference_normal_distributed('cylinderorientation', ['img', 'result', 'result', 1, 1], 'img', input_file(), convert_to_type=ImageDataType.FLOAT32, maxmem=100)
    
    def test_rot90cw_dist(self):
        self.check_difference_normal_distributed('rot90cw', ['img', 'result'], 'result', input_file(), maxmem=10)
    
    def test_rot90ccw_dist(self):
        self.check_difference_normal_distributed('rot90ccw', ['img', 'result'], 'result', input_file(), maxmem=10)
    
    def test_reslice1_dist(self):
        self.check_difference_normal_distributed('reslice', ['img', 'result', 'top'], 'result', input_file(), maxmem=3)
    
    def test_reslice2_dist(self):
        self.check_difference_normal_distributed('reslice', ['img', 'result', 'bottom'], 'result', input_file(), maxmem=3)
    
    def test_reslice3_dist(self):
        self.check_difference_normal_distributed('reslice', ['img', 'result', 'left'], 'result', input_file(), maxmem=3)
    
    def test_reslice4_dist(self):
        self.check_difference_normal_distributed('reslice', ['img', 'result', 'right'], 'result', input_file(), maxmem=3)
    
    def test_flip1_dist(self):
        self.check_difference_normal_distributed('flip', ['img', 0], 'img', input_file(), maxmem=5)
    
    def test_flip2_dist(self):
        self.check_difference_normal_distributed('flip', ['img', 1], 'img', input_file(), maxmem=5)
    
    def test_flip3_dist(self):
        self.check_difference_normal_distributed('flip', ['img', 2], 'img', input_file(), maxmem=5)


    ## General rotations involve interpolation so we expect images to match only to a small tolerance.
    def test_rotate1_dist(self):
        self.check_difference_normal_distributed('rotate', ['img', 'result', 30/180*3.14], 'result', input_file(), convert_to_type=ImageDataType.FLOAT32, tolerance=0.1, maxmem=10)
    
    def test_rotate2_dist(self):
        self.check_difference_normal_distributed('rotate', ['img', 'result', 30/180*3.14, [0, 0, 1], [128, 128, 64], [128, 128, 64]], 'result', input_file(), convert_to_type=ImageDataType.FLOAT32, tolerance=0.1, maxmem=10)
    
    def test_rotate3_dist(self):
        self.check_difference_normal_distributed('rotate', ['img', 'result', 30/180*3.14, [1, 0, 0], [128, 128, 64], [128, 128, 64]], 'result', input_file(), convert_to_type=ImageDataType.FLOAT32, tolerance=0.1, maxmem=20)
    
    def test_rotate4_dist(self):
        self.check_difference_normal_distributed('rotate', ['img', 'result', 30/180*3.14, [1, 1, 1], [128, 128, 64], [128, 128, 64]], 'result', input_file(), convert_to_type=ImageDataType.FLOAT32, tolerance=0.1, maxmem=40)
    
    def test_rotate5_dist(self):
        self.check_difference_normal_distributed('rotate', ['img', 'result', 30/180*3.14, [0, 0, 1], [10, 10, 0], [10, 10, 0]], 'result', input_file(), convert_to_type=ImageDataType.FLOAT32, tolerance=0.1, maxmem=10)
    
    def test_rotate6_dist(self):
        self.check_difference_normal_distributed('rotate', ['img', 'result', 30/180*3.14, [1, 0, 0], [10, 10, 0], [10, 10, 0]], 'result', input_file(), convert_to_type=ImageDataType.FLOAT32, tolerance=0.1, maxmem=40)
    
    def test_rotate7_dist(self):
        self.check_difference_normal_distributed('rotate', ['img', 'result', 30/180*3.14, [1, 1, 1], [10, 10, 0], [10, 10, 0]], 'result', input_file(), convert_to_type=ImageDataType.FLOAT32, tolerance=0.1, maxmem=40)
    
    def test_rotate8_dist(self):
        self.check_difference_normal_distributed('rotate', ['img', 'result', 30/180*3.14, [0, 0, 1], [10, 10, 0], [30, 20, 40]], 'result', input_file(), convert_to_type=ImageDataType.FLOAT32, tolerance=0.1, maxmem=10)
    
    def test_rotate9_dist(self):
        self.check_difference_normal_distributed('rotate', ['img', 'result', 30/180*3.14, [1, 0, 0], [10, 10, 0], [30, 20, 40]], 'result', input_file(), convert_to_type=ImageDataType.FLOAT32, tolerance=0.1, maxmem=20)
    
    def test_rotate10_dist(self):
        self.check_difference_normal_distributed('rotate', ['img', 'result', 30/180*3.14, [1, 1, 1], [10, 10, 0], [30, 20, 40]], 'result', input_file(), convert_to_type=ImageDataType.FLOAT32, tolerance=0.1, maxmem=40)
    
    def test_meancurvature_dist(self):
        self.check_difference_normal_distributed('meancurvature', ['img', 'result', 1], 'result', maxmem=50)

    ## Here we use tolerance as curvature seems to be sensitive to the origin of the calculation blocks in the distributed mode.
    def test_curvature_dist(self):
        self.check_difference_normal_distributed('curvature', ['img', 10, 'result', 'result'], 'result', input_file_bin(), convert_to_type=ImageDataType.FLOAT32, maxmem=30, tolerance=0.15)
    
    def test_set_dist(self):
        self.check_difference_normal_distributed('set', ['img', [100, 100, 100], 100], 'img', maxmem=5)
    
    def test_sphere_dist(self):
        self.check_difference_normal_distributed('sphere', ['img', [100, 100, 100], 100, 200], 'img', maxmem=5)
    
    def test_box_dist(self):
        self.check_difference_normal_distributed('box', ['img', [100, 100, 100], [100, 80, 60], 200], 'img', maxmem=5)
    
    def test_ellipsoid_dist(self):
        self.check_difference_normal_distributed('ellipsoid', ['img', [100, 100, 100], [100, 80, 60], 200, [1, 1, 1], [1, 1, -1]], 'img', maxmem=5)
    
    def test_box_dist(self):
        self.check_difference_normal_distributed('box', ['img', [100, 100, 100], [100, 80, 60], 200, [1, 1, 1], [1, 1, -1]], 'img', maxmem=5)
    
    def test_ramp1_dist(self):
        self.check_difference_normal_distributed('ramp', ['img', 2], 'img', maxmem=5)
    
    def test_ramp2_dist(self):
        self.check_difference_normal_distributed('ramp', ['img', 1], 'img', maxmem=5)
    
    def test_ramp3_dist(self):
        self.check_difference_normal_distributed('ramp', ['img', 0], 'img', maxmem=5)
    
    def test_line_dist(self):
        self.check_difference_normal_distributed('line', ['img', [10, 20, 30], [100, 200, 100], 200], 'img', maxmem=5)
    
    def test_capsule_dist(self):
        self.check_difference_normal_distributed('capsule', ['img', [10, 20, 30], [100, 200, 100], 10, 200], 'img', maxmem=5)
    
    def test_stddev_dist(self):
        self.check_difference_normal_distributed('stddev', ['img', 'result'], 'result', maxmem=5, tolerance=1e-4)
    
    def test_autothreshold_dist(self):
        self.check_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.OTSU], 'img', maxmem=5)
        self.check_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.HUANG], 'img', maxmem=5)
        self.check_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.ISODATA], 'img', maxmem=5)
        self.check_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.LI], 'img', maxmem=5)
        self.check_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.MAXENTROPY], 'img', maxmem=5)
        self.check_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.MEAN], 'img', maxmem=5)
        self.check_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.MINERROR], 'img', maxmem=5)
        self.check_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.MOMENTS], 'img', maxmem=5)
        self.check_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.PERCENTILE], 'img', maxmem=5)
        self.check_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.RENYI], 'img', maxmem=5)
        self.check_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.SHANBHAG], 'img', maxmem=5)
        self.check_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.TRIANGLE], 'img', maxmem=5)
        self.check_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.YEN], 'img', maxmem=5)
        self.check_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.MEDIAN], 'img', maxmem=5)
        self.check_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.MIDGREY], 'img', maxmem=5)
        self.check_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.NIBLACK], 'img', maxmem=5)
        self.check_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.PHANSALKAR], 'img', maxmem=5)
        self.check_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.SAUVOLA], 'img', maxmem=5)
        self.check_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.BERNSEN], 'img', maxmem=5)
    
    def test_autothresholdvalue_dist(self):
        self.check_difference_normal_distributed('autothresholdvalue', ['img', 'result', AutoThresholdMethod.OTSU], 'result', maxmem=5)
        self.check_difference_normal_distributed('autothresholdvalue', ['img', 'result', AutoThresholdMethod.HUANG], 'result', maxmem=5)
        self.check_difference_normal_distributed('autothresholdvalue', ['img', 'result', AutoThresholdMethod.ISODATA], 'result', maxmem=5)
        self.check_difference_normal_distributed('autothresholdvalue', ['img', 'result', AutoThresholdMethod.LI], 'result', maxmem=5)
        self.check_difference_normal_distributed('autothresholdvalue', ['img', 'result', AutoThresholdMethod.MAXENTROPY], 'result', maxmem=5)
        self.check_difference_normal_distributed('autothresholdvalue', ['img', 'result', AutoThresholdMethod.MEAN], 'result', maxmem=5)
        self.check_difference_normal_distributed('autothresholdvalue', ['img', 'result', AutoThresholdMethod.MINERROR], 'result', maxmem=5)
        self.check_difference_normal_distributed('autothresholdvalue', ['img', 'result', AutoThresholdMethod.MOMENTS], 'result', maxmem=5)
        self.check_difference_normal_distributed('autothresholdvalue', ['img', 'result', AutoThresholdMethod.PERCENTILE], 'result', maxmem=5)
        self.check_difference_normal_distributed('autothresholdvalue', ['img', 'result', AutoThresholdMethod.RENYI], 'result', maxmem=5)
        self.check_difference_normal_distributed('autothresholdvalue', ['img', 'result', AutoThresholdMethod.SHANBHAG], 'result', maxmem=5)
        self.check_difference_normal_distributed('autothresholdvalue', ['img', 'result', AutoThresholdMethod.TRIANGLE], 'result', maxmem=5)
        self.check_difference_normal_distributed('autothresholdvalue', ['img', 'result', AutoThresholdMethod.YEN], 'result', maxmem=5)
        self.check_difference_normal_distributed('autothresholdvalue', ['img', 'result', AutoThresholdMethod.MEDIAN], 'result', maxmem=5)
        self.check_difference_normal_distributed('autothresholdvalue', ['img', 'result', AutoThresholdMethod.MIDGREY], 'result', maxmem=5)
        self.check_difference_normal_distributed('autothresholdvalue', ['img', 'result', AutoThresholdMethod.NIBLACK], 'result', maxmem=5)
        self.check_difference_normal_distributed('autothresholdvalue', ['img', 'result', AutoThresholdMethod.PHANSALKAR], 'result', maxmem=5)
        self.check_difference_normal_distributed('autothresholdvalue', ['img', 'result', AutoThresholdMethod.SAUVOLA], 'result', maxmem=5)
        self.check_difference_normal_distributed('autothresholdvalue', ['img', 'result', AutoThresholdMethod.BERNSEN], 'result', maxmem=5)
        
    def test_autothresholdvalue_custom(self):
        pi = self.pi2
        img = pi.read(input_file())
        th = pi.autothresholdvalue(img)
        print(th)

    def test_sum_custom(self):
        pi = self.pi2
        img = pi.read(input_file())
        ret = pi.sum(img)
        print(ret)

    def test_mean_custom(self):
        pi = self.pi2
        img = pi.read(input_file())
        ret = pi.mean(img)
        print(ret)

    def test_stddev_custom(self):
        pi = self.pi2
        img = pi.read(input_file())
        ret = pi.stddev(img)
        print(ret)

    def test_minval_custom(self):
        pi = self.pi2
        img = pi.read(input_file())
        ret = pi.minval(img)
        print(ret)

    def test_maxval_custom(self):
        pi = self.pi2
        img = pi.read(input_file())
        ret = pi.maxval(img)
        print(ret)

    def test_localthreshold_dist(self):
        self.check_difference_normal_distributed('localthreshold', ['img', 'result', [2, 2, 2], AutoThresholdMethod.OTSU], 'result', maxmem=5)
    
    def test_add_dist(self):
        self.check_difference_normal_distributed('add', ['result', 'img'], 'result', maxmem=5)
    
    def test_eval1_dist(self):
        self.check_difference_normal_distributed('eval', ['77', 'img'], 'img', maxmem=5)
    
    def test_eval2_dist(self):
        self.check_difference_normal_distributed('eval', ['x0+x1', 'result', 'img'], 'result', maxmem=5)
    
    def test_eval3_dist(self):
        self.check_difference_normal_distributed('eval', ['x0+x1*x2', 'result', 'img', 'img'], 'result', maxmem=5)
    
    def test_ramp3_1_dist(self):
        self.check_difference_normal_distributed('ramp3', ['img'], 'img', maxmem=5)

    def test_isimagefile(self):
        pi = self.pi2
        self.check_result(pi.isimagefile(input_file()) == True, "input image is tested not to exist")
        self.check_result(pi.isimagefile('non-existing file name') == False, "not existing image is test to exist")
        pi.distribute(Distributor.Local)
        self.check_result(pi.isimagefile(input_file()) == True, "input image is tested not to exist (distributed)")
        self.check_result(pi.isimagefile('non-existing file name') == False, "not existing image is test to exist (distributed)")
        pi.distribute(Distributor.No);


    def test_delaying1(self):
        infile = input_file()
        self.check_difference_delaying('delaying_1', f"read(img, {infile}); gaussfilter(img, out, 1); add(img, 100); clear(img); convert(out, conv, uint8);", 'conv');

    def test_delaying1(self):
        infile = input_file()
        self.check_difference_delaying('delaying_2', f"read(img, {infile}); add(img, 100); crop(img, out, [0,0,0], [100, 100, 120]); subtract(out, 100);", 'out');
    
    def test_delaying1(self):
        infile = input_file()
        self.check_difference_delaying('delaying_3', f"read(img, {infile}); convert(img, img32, float32); clear(img); threshold(img32, 5e-4); convert(img32, cyl, uint8); clear(img32);", 'cyl');
    
    def test_delaying1(self):
        infile = input_file()
        self.check_difference_delaying('delaying_4', f"read(img, {infile}); convert(img, img32, float32); clear(img); cylindricality(img32, 0.5, 0.5); threshold(img32, 5e-4); convert(img32, cyl, uint8); clear(img32);", 'cyl', maxmem=100);

    def test_numpy(self):
        infile = input_file()
        pi = self.pi2

        img = pi.read(infile)
        print(img.dimensions())

        # This should work without exceptions about image dimensions.
        nparr = img.to_numpy()
        nparr += 1
        pi.add(img, nparr)


    def test_surfacearea(self):
        infile = input_file()
        pi = self.pi2

        img = pi.read(infile)
        area = pi.surfacearea(img, 160)
        print(f"Area = {area}")

        assert round(area) == 325995, "Incorrect surface area estimate."


    ## These commands are not guaranteed to give the same result in distributed and local processing
    ## mode so they are not tested here. Despite that, they should give a valid result in both cases, but their
    ## result is not unique.
    ##self.check_difference_normal_distributed('surfacethin', ['result'])
    ##self.check_difference_normal_distributed('surfaceskeleton', ['result'])
    ##self.check_difference_normal_distributed('linethin', ['result'])
    ##self.check_difference_normal_distributed('lineskeleton', ['result'])

    ## This test involves pretty large dataset, so it is skipped here.
    ##self.check_difference_normal_distributed('scalelabels', ['img', 'result', 4])

    ## This test will never succeed as the approximate bilateral filtering algorithm uses random numbers.
    ##self.check_difference_normal_distributed('bilateralfilterapprox', ['img', 'result', 5, 200], 'result') 

    ## These two tests do not succeed as the default settings assume the data is spread on the whole 16-bit value range
    ##self.check_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.INTERMODES], 'img', maxmem=5)
    ##self.check_difference_normal_distributed('autothreshold', ['img', AutoThresholdMethod.MINIMUM], 'img', maxmem=5)
