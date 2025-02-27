source 'https://rubygems.org'
ruby ">= 2.7", "< 4"

gem 'github-pages'
gem 'rouge'
gem "webrick", "~> 1.7"

# rdiscount is a dependency of github-pages but earlier
# version may trigger issue https://github.com/davidfstr/rdiscount/issues/145
# (at least on macOS) so we explicitly pin here
gem "rdiscount", ">= 2.2.0.2"

# Windows and JRuby does not include zoneinfo files, so bundle the tzinfo-data gem
# and associated library; see https://jekyllrb.com/docs/installation/windows/
platforms :mingw, :x64_mingw, :mswin, :jruby do
  gem "tzinfo", ">= 1", "< 3"
  gem "tzinfo-data"
end
