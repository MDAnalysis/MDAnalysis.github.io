# Running Jekyll locally

Follow the [GitHub: Testing your GitHub Pages site locally with
Jekyll](https://docs.github.com/en/pages/setting-up-a-github-pages-site-with-jekyll/testing-your-github-pages-site-locally-with-jekyll) instructions

Two possible approaches: Either install `ruby` (>= 2.7) and
[jekyll](https://jekyllrb.com/) locally or use a Docker image with a
jekyll installation. The local ruby/jekyll is more flexible if you can
get it to work.

## Local jekyll installation

Basically follow [GitHub: Testing your GitHub Pages site locally with
Jekyll](https://docs.github.com/en/pages/setting-up-a-github-pages-site-with-jekyll/testing-your-github-pages-site-locally-with-jekyll)
and in particular [Jekyll's  installation
instructions](https://jekyllrb.com/docs/installation/).


### macOS

Do not use the system ruby (i.e., not `/usr/bin/ruby`) but install
your own version. The [Jekyll macOS
instructions](https://jekyllrb.com/docs/installation/macos/) use
[homebrew](https://brew.sh/) but this works equally well with
[macports](https://www.macports.org/). 

#### Overview

1. Use `chruby` and `ruby-install` to manage different ruby versions.
2. Build ruby 3.4.1 (this was the recommended version in the Jekyll
   instructions at the time)
3. Install all necessary gems from [Gemfile][] with `bundle`.


#### ruby installation with macports

Install `chruby` and `ruby-install`: On macos with macports

    sudo port install chruby ruby-install

Update `~/.bash_profile` (or your `~/.zshrc`) with

    # chruby (macports)
    source /opt/local/share/chruby/chruby.sh

    # To enable auto-switching of Rubies specified by .ruby-version files, add the
    # following to ~/.bash_profile or ~/.zshrc:
    source /opt/local/share/chruby/auto.sh

Install ruby (compile from source, this requires Xcode and development
tools on macOS and will take a while):

    ruby-install ruby 3.4.1

Switch to using it (can also be added to `~/.bash_profile`)

    chruby ruby-3.4.1

To always activate this version, add a `.ruby-version` file in the
directory with the `Gemfile` (see [How to install different versions
of ruby and automatically switch between them](https://www.moncefbelyamani.com/how-to-install-xcode-homebrew-git-rvm-ruby-on-mac/#how-to-install-different-versions-of-ruby-and-switch-between-them)
-- also needs chruby to be active):

    echo "3.4.1" > .ruby-version

(May need to open a new terminal to let `.bash_profile` changes to take effect).


### All OS: install gems (including jekyll)

Ensure that the correct version of ruby is use (`ruby --version`) and
then 

    bundle install
	
should download and build the appropriate gems.

The `Gemfile.lock` should show something like
```ruby
GEM
  remote: https://rubygems.org/
  specs:
    ...
    github-pages (232)
      github-pages-health-check (= 1.18.2)
      jekyll (= 3.10.0)
      ...
```

### Run jekyll locally

Start jekyll and have it serve the site at http://127.0.0.1:4000
```bash
rm -r _site
bundle exec jekyll serve --watch --livereload --future --incremental
```

(See the [serve commandline
options](https://jekyllrb.com/docs/configuration/options/#serve-command-options)
for more details.)



## Docker image

Trying a stupid, very old docker image

### macOS docker installation

- install Docker Desktop https://docs.docker.com/desktop/install/mac-install/
- start docker
- follow the solution from  https://stackoverflow.com/a/58850151/ as
  described next


### Build the site

To build the static site, run `jekyll build` inside docker:
```bash
export JEKYLL_VERSION=3.8
docker run --rm \
  --volume="$PWD:/srv/jekyll" \
  -it jekyll/jekyll:$JEKYLL_VERSION \
  jekyll build
```
The static site files are stored in the `_site` directory.

### Serve the static site

Then you must *serve* the site:
```bash
cd _site
python -m http.serve
```

Point your browser to http://localhost:8000

When you make changes, you need to re-build the site.



