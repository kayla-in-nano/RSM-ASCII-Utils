"""
    ASCII-RSM-Functions.jl: defines various functions useful for extracting and plotting XRD RSM data from RINT ASCII raw datafiles.
    Copyright (C) 2025 Mikayla Lord (https://orcid.org/0000-0002-1388-3872)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    If you distribute this code or use it to process data, I would appreciate it if you attributed or cited me as the author.

    Contact me via email: kaylaTAL@protonmail.com
    Github: https://github.com/kayla-in-nano/RSM-ASCII-Utils

    If you distribute this code or 
"""

using DelimitedFiles
using DataFrames
using CSV
using LaTeXStrings
using PlotlyJS

"""
    qspace_filter(qx, qz, min_qx, max_qx, min_qz, max_qz)::Bool

Return the trimed dataframe that only includes rows min_qx < qx < max_qx and  min_qz < qz < max_qx.
"""
function qspace_filter(qx, qz, min_qx::Number, max_qx::Number, min_qz::Number, max_qz::Number)::Bool #this and gonio filter should be one single function that depends on the range of the final plot,
    interesting_qx = min_qx < qx && max_qx > qx # Condition for qx.  #and choosing to filter based on gonio or qspace coords should be decided by functions that use them 
    interesting_qz = min_qz < qz && max_qz > qz #Condition for qy.
    interesting_qx && interesting_qz #Trims dataframe so that only rows that satisfy the union of the conditions for qx and qy are included, may speed up plotting.
end

"""
    gonio_filter(x, twotheta, min_x, max_x, min_twotheta, max_twotheta)::Bool

Return the trimed dataframe that only includes rows where min_x < x < max_x and  min_twotheta < twotheta < max_twotheta.

Here, x is the other axis that should be plotted on the x axis, such as omega or chi.
"""
function gonio_filter(x, twotheta, min_x::Number, max_x::Number, min_twotheta::Number, max_twotheta::Number)::Bool
    interesting_x = min_x < x && max_x > x
    interesting_2theta = min_twotheta < twotheta && max_twotheta > twotheta
    interesting_x && interesting_2theta
end


"""
    gen_dir_string(hkl::String; sub = nothing)

Format using markdown a string of three integers separated by spaces as a crystallographic direction.

The hkl string should include only integer numbers for each direction (positive or negative), and spaces to separate them. Negative indicies are given overlines.

It is not necessary to include brackets. An optional sub argument can be specified as a string, and will be included as a subscript to the direction.

# Examples
```jldoctest
julia> gen_dir_string("1 -2 3", sub="pc")
"[1<span style='text-decoration:overline'>2</span>3]<sub>pc</sub>"
```
"""
function gen_dir_string(hkl::String; sub = nothing)
    hkl_int = parse.(Int, split(hkl, ' ')) #parse the hkl string as a list delimited by spaces (this function could've been written to take a list as input but I found that clunky when in use)
    hkl_str = Vector{String}() #List of strings to hold h, k, and l after formatting.
    for i in hkl_int
        if i < 0 #case for negative index
            i = i * -1 #make number positive
            push!(hkl_str, "<span style='text-decoration:overline'>$i</span>") #add an overline, then add to hkl_str
        else #case for positive number
            push!(hkl_str, "$i") #Convert Int to String, then add to hkl_str
        end
    end
    if isnothing(sub) == true #no subscript
        dir_str = "["*hkl_str[1]*hkl_str[2]*hkl_str[3]*"]" #Combine formatted h, k, and l, and add brackets
    else #there is a subscript
        dir_str = "["*hkl_str[1]*hkl_str[2]*hkl_str[3]*"]<sub>$sub</sub>" #Combine formatted h, k, and l, add brackets, and add the subscript
    end
    return dir_str
end

"""
    gen_Qaxis_title(hkl::String; sub = nothing)

Format using markdown an axis title of the form "Q \\\\ [h k l]_sub (Å^-1)".

The hkl string should include only integer numbers for each direction (positive or negative), and spaces to separate them. Negative indicies are given overlines.

It is not necessary to include brackets. An optional sub argument can be specified as a string, and will be included as a subscript to the direction.

See also [`gen_dir_string`](@ref).

# Examples
```jldoctest
julia> gen_Qaxis_title("1 -2 3", sub="pc")
"<i><b>Q</b></i>&#8201;//&#8201;[1<span style='text-decoration:overline'>2</span>3]<sub>pc</sub> (Å<sup>-1</sup>)"
```
"""
function gen_Qaxis_title(hkl::String; sub = nothing)
    if isnothing(sub) == true 
        direction = gen_dir_string(hkl)
    else
        direction = gen_dir_string(hkl, sub = sub)
    end
    axis_title = "<i><b>Q</b></i>&#8201;//&#8201;"*direction*" (Å<sup>-1</sup>)" #&#8201 is a Thin Space character
    return axis_title
end

#"Turbo" colourmap that transitions to white at low values. Default colourmap for plotting.
cmap = [[0.0, "rgb(255, 255, 255)"], [0.003937007874015748, "rgb(191, 183, 198)"], [0.007874015748031496, "rgb(128, 112, 142)"], [0.011811023622047244, "rgb(65, 40, 86)"], [0.015748031496062992, "rgb(51, 26, 79)"], [0.01968503937007874, "rgb(52, 29, 85)"], [0.023622047244094488, "rgb(53, 31, 92)"], [0.027559055118110236, "rgb(54, 34, 99)"], [0.031496062992125984, "rgb(56, 38, 107)"], [0.03543307086614173, "rgb(57, 42, 115)"], [0.03937007874015748, "rgb(58, 45, 123)"], [0.04330708661417323, "rgb(59, 48, 128)"], [0.047244094488188976, "rgb(60, 50, 134)"], [0.051181102362204724, "rgb(61, 53, 139)"], [0.05511811023622047, "rgb(61, 55, 145)"], [0.05905511811023622, "rgb(62, 58, 150)"], [0.06299212598425197, "rgb(63, 61, 155)"], [0.06692913385826771, "rgb(63, 63, 160)"], [0.07086614173228346, "rgb(64, 66, 165)"], [0.07480314960629922, "rgb(65, 68, 170)"], [0.07874015748031496, "rgb(65, 72, 175)"], [0.0826771653543307, "rgb(66, 75, 181)"], [0.08661417322834646, "rgb(67, 78, 187)"], [0.09055118110236221, "rgb(67, 81, 192)"], [0.09448818897637795, "rgb(68, 84, 195)"], [0.09842519685039369, "rgb(68, 86, 199)"], [0.10236220472440945, "rgb(68, 89, 203)"], [0.1062992125984252, "rgb(69, 91, 207)"], [0.11023622047244094, "rgb(69, 94, 210)"], [0.11417322834645668, "rgb(69, 96, 214)"], [0.11811023622047244, "rgb(69, 99, 217)"], [0.1220472440944882, "rgb(70, 101, 220)"], [0.12598425196850394, "rgb(70, 104, 223)"], [0.12992125984251968, "rgb(70, 107, 226)"], [0.13385826771653542, "rgb(70, 110, 230)"], [0.1377952755905512, "rgb(70, 113, 233)"], [0.14173228346456693, "rgb(70, 116, 236)"], [0.14566929133858267, "rgb(70, 118, 238)"], [0.14960629921259844, "rgb(70, 121, 240)"], [0.15354330708661418, "rgb(70, 123, 242)"], [0.15748031496062992, "rgb(70, 125, 244)"], [0.16141732283464566, "rgb(70, 128, 246)"], [0.1653543307086614, "rgb(70, 130, 248)"], [0.16929133858267717, "rgb(69, 132, 249)"], [0.1732283464566929, "rgb(69, 134, 250)"], [0.17716535433070865, "rgb(69, 137, 252)"], [0.18110236220472442, "rgb(68, 140, 252)"], [0.18503937007874016, "rgb(67, 143, 253)"], [0.1889763779527559, "rgb(66, 146, 254)"], [0.19291338582677164, "rgb(65, 149, 254)"], [0.19685039370078738, "rgb(64, 151, 254)"], [0.20078740157480315, "rgb(63, 153, 254)"], [0.2047244094488189, "rgb(61, 156, 253)"], [0.20866141732283464, "rgb(60, 158, 253)"], [0.2125984251968504, "rgb(59, 160, 252)"], [0.21653543307086615, "rgb(57, 163, 251)"], [0.2204724409448819, "rgb(56, 165, 250)"], [0.22440944881889763, "rgb(54, 167, 249)"], [0.22834645669291337, "rgb(52, 170, 248)"], [0.23228346456692914, "rgb(50, 173, 246)"], [0.23622047244094488, "rgb(48, 176, 244)"], [0.24015748031496062, "rgb(46, 179, 242)"], [0.2440944881889764, "rgb(44, 181, 240)"], [0.24803149606299213, "rgb(42, 184, 238)"], [0.25196850393700787, "rgb(41, 186, 236)"], [0.2559055118110236, "rgb(39, 188, 234)"], [0.25984251968503935, "rgb(38, 190, 232)"], [0.2637795275590551, "rgb(36, 192, 230)"], [0.26771653543307083, "rgb(35, 194, 227)"], [0.27165354330708663, "rgb(33, 197, 225)"], [0.2755905511811024, "rgb(32, 199, 223)"], [0.2795275590551181, "rgb(31, 201, 220)"], [0.28346456692913385, "rgb(29, 203, 217)"], [0.2874015748031496, "rgb(28, 206, 214)"], [0.29133858267716534, "rgb(26, 208, 211)"], [0.2952755905511811, "rgb(25, 210, 208)"], [0.2992125984251969, "rgb(25, 212, 205)"], [0.3031496062992126, "rgb(24, 214, 203)"], [0.30708661417322836, "rgb(24, 216, 200)"], [0.3110236220472441, "rgb(23, 217, 198)"], [0.31496062992125984, "rgb(23, 219, 196)"], [0.3188976377952756, "rgb(23, 221, 193)"], [0.3228346456692913, "rgb(24, 222, 191)"], [0.32677165354330706, "rgb(24, 224, 189)"], [0.3307086614173228, "rgb(24, 225, 186)"], [0.3346456692913386, "rgb(26, 227, 183)"], [0.33858267716535434, "rgb(27, 229, 180)"], [0.3425196850393701, "rgb(28, 230, 178)"], [0.3464566929133858, "rgb(30, 232, 175)"], [0.35039370078740156, "rgb(32, 233, 173)"], [0.3543307086614173, "rgb(33, 234, 170)"], [0.35826771653543305, "rgb(35, 235, 168)"], [0.36220472440944884, "rgb(38, 236, 165)"], [0.3661417322834646, "rgb(40, 238, 162)"], [0.3700787401574803, "rgb(42, 239, 159)"], [0.37401574803149606, "rgb(45, 240, 156)"], [0.3779527559055118, "rgb(47, 241, 154)"], [0.38188976377952755, "rgb(50, 242, 150)"], [0.3858267716535433, "rgb(54, 243, 146)"], [0.38976377952755903, "rgb(58, 244, 142)"], [0.39370078740157477, "rgb(62, 245, 138)"], [0.39763779527559057, "rgb(66, 246, 135)"], [0.4015748031496063, "rgb(69, 247, 132)"], [0.40551181102362205, "rgb(73, 248, 129)"], [0.4094488188976378, "rgb(76, 248, 126)"], [0.41338582677165353, "rgb(80, 249, 122)"], [0.41732283464566927, "rgb(83, 250, 119)"], [0.421259842519685, "rgb(87, 250, 116)"], [0.4251968503937008, "rgb(90, 251, 113)"], [0.42913385826771655, "rgb(94, 251, 110)"], [0.4330708661417323, "rgb(98, 252, 107)"], [0.43700787401574803, "rgb(103, 252, 103)"], [0.4409448818897638, "rgb(108, 253, 99)"], [0.4448818897637795, "rgb(113, 253, 95)"], [0.44881889763779526, "rgb(116, 254, 92)"], [0.452755905511811, "rgb(120, 254, 89)"], [0.45669291338582674, "rgb(123, 254, 86)"], [0.46062992125984253, "rgb(127, 254, 84)"], [0.4645669291338583, "rgb(130, 254, 81)"], [0.468503937007874, "rgb(134, 254, 79)"], [0.47244094488188976, "rgb(137, 254, 76)"], [0.4763779527559055, "rgb(140, 254, 74)"], [0.48031496062992124, "rgb(144, 254, 72)"], [0.484251968503937, "rgb(147, 254, 69)"], [0.4881889763779528, "rgb(151, 254, 67)"], [0.4921259842519685, "rgb(155, 253, 64)"], [0.49606299212598426, "rgb(158, 253, 62)"], [0.5, "rgb(163, 252, 59)"], [0.5039370078740157, "rgb(166, 251, 58)"], [0.5078740157480315, "rgb(168, 251, 57)"], [0.5118110236220472, "rgb(171, 250, 56)"], [0.515748031496063, "rgb(173, 249, 55)"], [0.5196850393700787, "rgb(176, 249, 54)"], [0.5236220472440944, "rgb(178, 248, 54)"], [0.5275590551181102, "rgb(181, 247, 53)"], [0.5314960629921259, "rgb(183, 246, 52)"], [0.5354330708661417, "rgb(187, 245, 52)"], [0.5393700787401575, "rgb(190, 243, 52)"], [0.5433070866141733, "rgb(193, 242, 51)"], [0.547244094488189, "rgb(195, 241, 51)"], [0.5511811023622047, "rgb(197, 239, 51)"], [0.5551181102362205, "rgb(200, 238, 51)"], [0.5590551181102362, "rgb(202, 237, 52)"], [0.562992125984252, "rgb(204, 235, 52)"], [0.5669291338582677, "rgb(207, 234, 52)"], [0.5708661417322834, "rgb(209, 233, 52)"], [0.5748031496062992, "rgb(211, 231, 52)"], [0.5787401574803149, "rgb(213, 230, 53)"], [0.5826771653543307, "rgb(215, 228, 53)"], [0.5866141732283464, "rgb(218, 226, 54)"], [0.5905511811023622, "rgb(221, 224, 54)"], [0.5944881889763779, "rgb(223, 222, 55)"], [0.5984251968503937, "rgb(225, 220, 55)"], [0.6023622047244095, "rgb(227, 218, 55)"], [0.6062992125984252, "rgb(229, 216, 56)"], [0.610236220472441, "rgb(230, 215, 56)"], [0.6141732283464567, "rgb(232, 213, 56)"], [0.6181102362204725, "rgb(234, 211, 57)"], [0.6220472440944882, "rgb(235, 209, 57)"], [0.6259842519685039, "rgb(237, 207, 57)"], [0.6299212598425197, "rgb(238, 205, 57)"], [0.6338582677165354, "rgb(240, 203, 58)"], [0.6377952755905512, "rgb(242, 200, 58)"], [0.6417322834645669, "rgb(243, 198, 58)"], [0.6456692913385826, "rgb(245, 196, 58)"], [0.6496062992125984, "rgb(246, 194, 57)"], [0.6535433070866141, "rgb(247, 192, 57)"], [0.6574803149606299, "rgb(248, 190, 57)"], [0.6614173228346456, "rgb(249, 188, 56)"], [0.6653543307086615, "rgb(249, 186, 56)"], [0.6692913385826772, "rgb(250, 184, 55)"], [0.6732283464566929, "rgb(251, 181, 55)"], [0.6771653543307087, "rgb(251, 179, 54)"], [0.6811023622047244, "rgb(252, 177, 53)"], [0.6850393700787402, "rgb(252, 174, 52)"], [0.6889763779527559, "rgb(253, 171, 51)"], [0.6929133858267716, "rgb(253, 168, 50)"], [0.6968503937007874, "rgb(253, 165, 49)"], [0.7007874015748031, "rgb(253, 162, 48)"], [0.7047244094488189, "rgb(254, 160, 47)"], [0.7086614173228346, "rgb(254, 157, 46)"], [0.7125984251968503, "rgb(254, 155, 45)"], [0.7165354330708661, "rgb(254, 152, 44)"], [0.7204724409448818, "rgb(253, 149, 43)"], [0.7244094488188977, "rgb(253, 146, 41)"], [0.7283464566929134, "rgb(253, 144, 40)"], [0.7322834645669292, "rgb(253, 141, 39)"], [0.7362204724409449, "rgb(252, 137, 38)"], [0.7401574803149606, "rgb(252, 134, 36)"], [0.7440944881889764, "rgb(251, 130, 34)"], [0.7480314960629921, "rgb(251, 127, 33)"], [0.7519685039370079, "rgb(250, 124, 32)"], [0.7559055118110236, "rgb(249, 121, 30)"], [0.7598425196850394, "rgb(249, 118, 29)"], [0.7637795275590551, "rgb(248, 115, 28)"], [0.7677165354330708, "rgb(247, 113, 27)"], [0.7716535433070866, "rgb(247, 110, 26)"], [0.7755905511811023, "rgb(246, 107, 24)"], [0.7795275590551181, "rgb(245, 104, 23)"], [0.7834645669291338, "rgb(244, 102, 22)"], [0.7874015748031495, "rgb(243, 98, 21)"], [0.7913385826771654, "rgb(241, 95, 19)"], [0.7952755905511811, "rgb(240, 91, 18)"], [0.7992125984251969, "rgb(239, 89, 17)"], [0.8031496062992126, "rgb(238, 86, 16)"], [0.8070866141732284, "rgb(236, 84, 15)"], [0.8110236220472441, "rgb(235, 81, 14)"], [0.8149606299212598, "rgb(234, 79, 13)"], [0.8188976377952756, "rgb(233, 77, 12)"], [0.8228346456692913, "rgb(231, 75, 12)"], [0.8267716535433071, "rgb(230, 73, 11)"], [0.8307086614173228, "rgb(229, 70, 10)"], [0.8346456692913385, "rgb(227, 68, 10)"], [0.8385826771653543, "rgb(225, 66, 9)"], [0.84251968503937, "rgb(223, 63, 8)"], [0.8464566929133858, "rgb(221, 61, 8)"], [0.8503937007874016, "rgb(220, 59, 7)"], [0.8543307086614174, "rgb(218, 57, 7)"], [0.8582677165354331, "rgb(216, 55, 6)"], [0.8622047244094488, "rgb(214, 53, 6)"], [0.8661417322834646, "rgb(213, 51, 5)"], [0.8700787401574803, "rgb(211, 50, 5)"], [0.8740157480314961, "rgb(209, 48, 4)"], [0.8779527559055118, "rgb(207, 46, 4)"], [0.8818897637795275, "rgb(205, 45, 4)"], [0.8858267716535433, "rgb(203, 43, 3)"], [0.889763779527559, "rgb(201, 41, 3)"], [0.8937007874015748, "rgb(198, 39, 3)"], [0.8976377952755905, "rgb(195, 37, 2)"], [0.9015748031496063, "rgb(193, 35, 2)"], [0.905511811023622, "rgb(190, 33, 2)"], [0.9094488188976377, "rgb(188, 32, 2)"], [0.9133858267716535, "rgb(186, 30, 1)"], [0.9173228346456693, "rgb(184, 29, 1)"], [0.9212598425196851, "rgb(181, 28, 1)"], [0.9251968503937008, "rgb(179, 26, 1)"], [0.9291338582677166, "rgb(176, 25, 1)"], [0.9330708661417323, "rgb(174, 24, 1)"], [0.937007874015748, "rgb(171, 22, 1)"], [0.9409448818897638, "rgb(168, 20, 1)"], [0.9448818897637795, "rgb(164, 19, 1)"], [0.9488188976377953, "rgb(161, 17, 1)"], [0.952755905511811, "rgb(158, 16, 1)"], [0.9566929133858267, "rgb(155, 15, 1)"], [0.9606299212598425, "rgb(152, 14, 1)"], [0.9645669291338582, "rgb(149, 13, 1)"], [0.968503937007874, "rgb(146, 11, 1)"], [0.9724409448818897, "rgb(144, 10, 1)"], [0.9763779527559056, "rgb(141, 9, 1)"], [0.9803149606299213, "rgb(137, 8, 1)"], [0.984251968503937, "rgb(134, 7, 1)"], [0.9881889763779528, "rgb(131, 6, 2)"], [0.9921259842519685, "rgb(128, 5, 2)"], [0.9960629921259843, "rgb(125, 4, 2)"], [1.0, "rgb(122, 4, 2)"]]



"""
    plot_gonio_map(ascifile="string"; cmap=cmap, xrange = nothing, yrange = nothing, offset_omega::Float64 = 0.0, offset_2theta::Float64 = 0.0)

Given a path to an RSM ASCII raw data file, returns a PlotlyJS plot, mapping the RSM in angular coordinates, and a one row dataframe giving details of the datapoint with the most counts.

Currently supports scans of type "2theta" (vs omega) and "2Theta/Omega" (vs omega). Support for 2Theta/Omega vs Chi is planned.

# Arguments
- `xrange`: not implemented.
- `xrange`: not implemented.
- `offset_omega::Number = 0.0`: Offset to be subtracted from omega if correction is needed.
- `offset_2theta::Number = 0.0`: Offset to be subtracted from 2Theta if correction is needed.

"""
function plot_gonio_map(ascifile="string"; cmap=cmap, xrange = nothing, yrange = nothing, offset_omega::Number = 0.0, offset_2theta::Number = 0.0)
    regexSTART = r"(?:\*START\t\t=  )([^\r][^\n]+)(?:\n\*STOP\t\t=  )([^\n]+)(?:\n\*STEP\t\t=  )([^\n]+)(?:\n\*OFFSET\t\t=  )([^\n]+)"
    f = open(ascifile)

    filestring = read(f, String) #convert the full file to a string

    m = match(regexSTART, filestring) #apply regex to the file string

    start_theta = match(r"([\d.]+)",m.captures[1]).captures[1]
    end_theta = match(r"([\d.]+)",m.captures[2]).captures[1]
    step_theta = match(r"([\d.]+)",m.captures[3]).captures[1]
    
    
    scan_axis = match(r"(?:\*SCAN_AXIS\t=  )([^\r][^\n\r]+)", filestring).captures[1]
    print("Scan axis = " * scan_axis * " ")

    #Print the scan range and step size
    print("start = " * start_theta * " ")
    print("end = " * end_theta * " ")
    print("step size = " * step_theta * " ")
    print("First offset = " * match(r"([\d.]+)",m.captures[4]).captures[1] * " ") #This is only the first offset

    #Convert the scan limits and step size to floats instead of strings
    start_theta = parse.(Float64, start_theta)
    end_theta = parse.(Float64, end_theta)
    step_theta = parse.(Float64, step_theta)

    close(f) #close the version of the file that was a string

    #Find data

    #regexDATA = r"(^[\d][\d., \n]+)"m
    regexDATA = r"(?:\*COUNT\t\t=  )([^\r][^\n]+)(?:\n)([\d][\d., \n\re-]+)"s #Regex to find datapoints
    matches = eachmatch(regexDATA, filestring)
    matchcount = 0
    ndatapoints = Vector{Int64}()
    data_string_array = Vector{String}()
    for match in matches
        matchcount += 1
        push!(ndatapoints, parse.(Int64, match.captures[1])) #create array of the number of points in each scan
        push!(data_string_array, match.captures[2]) #create an array where each element is a string of each scan
    end
    #print(matchcount) #should be number of scans
    #print(ndatapoints[1]) # should be number of points in each scan
    points = LinRange(start_theta, end_theta, ndatapoints[1]) #generate the values of 2theta (assuming all scans have the same range)
    scan_df = DataFrame(scan = points) #Convert this into a dataframe. Name the first coloumn scan to get good indexing of the column names for later, this gets changed to :TwoTheta later on

    regexScan = r"([\d\., e-]+)" #separate lines of data in each scan

    for i in 1:matchcount
        scan_string_array = "" #empty string to add to
        linematches = eachmatch(regexScan, data_string_array[i])
        for match in linematches
            scan_string_array = scan_string_array * match.captures[1] * ", " #create one long string of the whole scan
        end
        scan_string_array = rstrip(scan_string_array, [' ',',']) #remove trailing delimiter

        #Treat the scan as a CSV string and add it as a column to the dataframe
        insertcols!(scan_df, i + 1, :scan => CSV.File(IOBuffer(scan_string_array); header = false, transpose = true).Column1; makeunique = true) #this is giving an error
    end

    rename!(scan_df, Dict(:scan => "TwoTheta")) #Change first column to :TwoTheta

    f = open(ascifile) #reopen new instance of the file
    stringarray = collect(eachline(f)) #Create an array of strings, where each line of the file is a new element in the array
    close(f) #close the file

    #print(ndims(stringarray))
    #print(". Length of the string array is ")
    #print(length(stringarray))

    offsetregex = r"(?:\*OFFSET\t\t=  )([^\r][^\n]+)" #regex to find offsets associated with each scan
    offsetarray = Vector{Float64}() #empty array to store offsets
    for s in stringarray
        m = match(offsetregex, s)
        if typeof(m) == RegexMatch #ignore any lines where there is no offset
            push!(offsetarray, parse.(Float64, m.captures[1])) #add the offset to the array as a float
        end
    end

    omega_df = DataFrame() # create empty dataframe for omega values

    #scans_omega = names(scan_df, Not(:TwoTheta)) #get a string of scan names

    scans_omega = Vector{Symbol}()
    for i in 1:length(names(scan_df, Not(:TwoTheta)))
        push!(scans_omega, Symbol(names(scan_df, Not(:TwoTheta))[i])) #copy scan names to the omega dataframe from the dataframe with number of counts
    end

    omega_df.TwoTheta = scan_df.TwoTheta #make the first column of the omega dataframe 2theta
    
    if scan_axis == "2theta" #when doing 2theta vs omega to create an RSM, offset values in the ascii are equal to omega.
        for i in eachindex(offsetarray)
            omega_df[:, scans_omega[i]] .= offsetarray[i]
        end
    else
        for i in eachindex(offsetarray) #when doing coupled 2theta/omega vs omega to create an RSM, the offset given by the ascii + theta = omega
            transform!(omega_df, [:TwoTheta] => ByRow((x) -> ((x * 0.5) + offsetarray[i])) => scans_omega[i])
        end
    end

    long_df = stack(omega_df, Not(:TwoTheta)) #reformat the dataframe to have omega as a column
    rename!(long_df, Dict(:variable => "Scan_No", :value => "Omega")) #Change first column to :TwoTheta

    long_df.Counts = stack(scan_df, Not(:TwoTheta)).value #add the counts for each point in the RSM
    transform!(long_df, [:Omega] => ByRow((x) -> (x-offset_omega)) => :Omega) #Correct omega by using a user defined omega offset
    transform!(long_df, [:TwoTheta] => ByRow((x) -> (x-offset_2theta)) => :TwoTheta) #Correct 2theta by using a user defined 2theta offset
    #To do: add option for counts per second
    transform!(long_df, [:Counts] => ByRow((x) -> (log10(x+0.0001))) => :LogCounts) #Create a column that is the log of the counts, to allow for display using a log scale
    transform!(long_df, [:Counts] => ByRow((x) -> (string(round(x, sigdigits=3)))) => :z_formatted) #create a coloumn with rounded z values so that we can hover over points to see the number of counts

    long_df.Countslog = replace(long_df.LogCounts, -Inf => -0.3) #replace missing data points with a really small number to avoid issues when plotting

    maxcounts = long_df[long_df.Counts .== maximum(long_df.Counts), :] #details of the point with the most counts, useful for determining if an offset needs to be applied
    print(maxcounts) #print these details (this should be optional)

    hovertemplate = "<b>ω: %{x:.3f} <br>2θ: %{y:.3f}</b><BR>Counts: %{customdata}<br>" #what should be displaued when we hover over points in the graph?
    zmin = 0.5 #colourbar min
    zmax = 4 #colourbar max. These should be arguments.
    p = plot(scattergl( #plotting is done as a scatter plot because it is faster and better represents the actual data points that were collected without interpolation (unlike a contour plot)
                        #during testing, plotting as a heatmap or contour plot would inflate the number of counts compared to a scatter plot, and imply greater resolution that is actually present.
        long_df, x=:Omega, y=:TwoTheta, customdata = :z_formatted, #customdata is fed to the mouse hover diplay
            marker=attr(color=:LogCounts, colorscale=cmap, showscale=true, symbol="diamond", size=1.8, cmin = zmin, cmax = zmax, #size may need to be tweaked depending on the range
                colorbar=attr(
                    title_text="Counts<br><sup>&nbsp;</sup>", #empty subscript to add padding between the title and the colourbar
                    ticks="outside",
                    tickmode= "array",
                    tickvals = collect(ceil(zmin):floor(zmax)),
                    ticktext = ["10", "10<sup>2</sup>", "10<sup>3</sup>", "10<sup>4</sup>",], #this will not work with different z ranges, needs to be fixed

                    )),

        mode="markers",
        hovertemplate=hovertemplate, #ticks = "outside",

    ),config=PlotConfig(scrollZoom=false, toImageButtonOptions=attr(
            format="png", # one of png, svg, jpeg, webp
            filename=splitext(basename(ascifile))[1],
            #height=900,
            #width=1400,
            scale=2 # Multiply title/legend/axis/canvas sizes by this factor
        ).fields))

    relayout!(
        p, 
        autosize=0,
        width=width+11,#900 #need to figure out a more reliable to make te plot square
        height=floor(Int,width * l_0 / l_1), #800
        #xaxis_range=[minx, maxx],
        #yaxis_range=[minz, maxz], #auto determine ranges when plotting in gonio coords
        plot_bgcolor="White",
        yaxis=attr(title_text="2θ", tickformat=".2f", showline=true, linewidth=1, linecolor="black", mirror=true, ticks = "outside"),
        xaxis=attr(title_text="ω", tickformat=".2f", showline=true, linewidth=1, linecolor="black", mirror=true, ticks = "outside"), #add option for other axes, such as chi
        font=attr(
            family="Arial",
            size=21,
            color="Black"
        ),
        coloraxis_colorbar=attr(title_text="test"),
        margin=attr(
            l=50,
            r=50,
            b=100,
            t=100,
            pad=2
        ),
        
    )
    return p, maxcounts #add option to only return p
end

"""
    compare_plot_to_expected_angles(plot, maxcounts, calc_2theta::Number, calc_omega::Number; show_plot::Bool = true, return_offsets::Bool = true)

Given a plotly plot, maxcounts datafile, and expected or calculated values of 2theta and omega, prints the difference between the expected coordinates and the coordinates of the point with most counts.

If show_plot is true, it will return a plot with lines indicating the expected coordinates. If return_offsets is true, returns the differences between the expected and actual coordinates.
"""
function compare_plot_to_expected_angles(plot, maxcounts, calc_2theta::Number, calc_omega::Number; show_plot::Bool = true, return_offsets::Bool = true)

    max_omega = maxcounts[1, :Omega]
    max_2theta = maxcounts[1, :TwoTheta]

    omega_offset = max_omega - calc_omega
    twotheta_offset =  max_2theta - calc_2theta

    print("Point with maximum counts has a difference between expected omega and omega of ", omega_offset)
    print("\nand between expected 2theta and 2theta of", twotheta_offset, "\n")

    if show_plot == true
        add_hline!(plot,calc_2theta)
        add_vline!(plot,calc_omega)
    end

    if show_plot == true && return_offsets == true
        return plot, omega_offset, twotheta_offset
    elseif show_plot == false && return_offsets == true
        return omega_offset, twotheta_offset
    elseif show_plot == true && return_offsets == false
        return plot
    end
end



"""
    compare_raw_to_expected(ascifile, calc_2theta::Number, calc_omega::Number; cmap=cmap, apply_offsets::Bool = false)

Given a raw ascii datafile and expected angular coordinates of the point with maximum counts, returns a plot of that datafile with the expected coordinates marked, and the difference between the expected and actual coordinates.

Output is a tuple in the form `p, omega_offset::Float64, twotheta_offset::Float64` where p is a PlotlyJS plot. If apply_offsets is true, then the plot will show the data after each offset is subtracted. Colormap of the plot can be set using cmap.

See also [`compare_plot_to_expected_angles`](@ref) and [`plot_corrected_RSM_ascii`](@ref).
"""
function compare_raw_to_expected(ascifile, calc_2theta::Number, calc_omega::Number; cmap=cmap, apply_offsets::Bool = false)
    p, maxcounts = plot_corrected_RSM_ascii(ascifile, cmap=cmap)
    if apply_offsets == true
        omega_offset, twotheta_offset = compare_plot_to_expected_angles(p, maxcounts, calc_2theta, calc_omega, show_plot = false, return_offsets = true)
        p = plot_corrected_RSM_ascii(ascifile, cmap=cmap, offset_omega = omega_offset, offset_2theta = twotheta_offset)
        add_hline!(p[1],calc_2theta)
        add_vline!(p[1],calc_omega)
        return p[1], omega_offset, twotheta_offset
    else
        p, omega_offset, twotheta_offset = compare_plot_to_expected_angles(p, maxcounts, calc_2theta, calc_omega, show_plot = true, return_offsets = true)
        return p, omega_offset, twotheta_offset
    end
end