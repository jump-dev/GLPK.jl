using Test
import GLPK

function glpk_tst_5()
    prev_term_out = GLPK.term_out(GLPK.OFF)

    datadir = joinpath(Pkg.dir(), "GLPK", "test", "data")

    dataf1 = joinpath(datadir, "data1.txt")
    dataf2 = joinpath(datadir, "data2.txt")

    data = GLPK.sdf_open_file(dataf1)
    sum = 0
    for i = 1:10
        num = GLPK.sdf_read_int(data)
        sum += num
    end
    @test sum == 3566
    GLPK.sdf_close_file(data)

    data = GLPK.sdf_open_file(dataf2)

    GLPK.sdf_read_int(data)
    GLPK.sdf_line(data)
    GLPK.sdf_read_num(data)
    GLPK.sdf_line(data)
    GLPK.sdf_read_item(data)
    GLPK.sdf_line(data)
    GLPK.sdf_read_item(data)
    GLPK.sdf_line(data)
    GLPK.sdf_read_text(data)
    GLPK.sdf_line(data)
    GLPK.sdf_read_text(data)
    GLPK.sdf_line(data)
    GLPK.sdf_read_text(data)
    GLPK.sdf_line(data)

    GLPK.sdf_close_file(data)

    GLPK.term_out(prev_term_out)
end

glpk_tst_5()
