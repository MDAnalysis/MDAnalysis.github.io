source 'https://rubygems.org'
ruby "~> 2.6"

gem 'github-pages'
gem 'rouge'

# https://github.com/github/pages-gem/issues/887#issuecomment-1773067676
install_if -> { ENV["GITHUB_ACTIONS"] != "true" } do
    puts "Is GitHub action: #{ENV["GITHUB_ACTIONS"] == "true"}"
    gem "webrick", "~> 1.7"
end 
