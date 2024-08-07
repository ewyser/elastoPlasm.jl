@info "generate README.md file.."
# Define the files to be concatenated
files = [
    "./docs/src/index.md",
]
# Open the README.md file for writing
open("README.md", "w") do readme
    for file ∈ files
        # Write the content of each file to README.md
        open(file, "r") do f
            content = read(f, String)
            write(readme, content)
            write(readme, "\n\n") # Add some space between sections
        end
    end
end
println("README.md has been generated successfully.")