import { NgModule } from '@angular/core';
import { BrowserModule } from '@angular/platform-browser';
import { AppRoutingModule } from './app-routing.module';
import { AppComponent } from './app.component';
import { FormsModule, ReactiveFormsModule } from '@angular/forms';
import { HttpClientModule } from '@angular/common/http';
import {NgxPaginationModule} from 'ngx-pagination'; 

// General imports
import { HomeComponent } from './layout/home/home.component';
import { FooterComponent } from './layout/footer/footer.component';
import { NavbarComponent } from './layout/navbar/navbar.component';
import { HeaderComponent } from './layout/header/header.component';
import { NotFoundComponent } from './meta/not-found/not-found.component';
import {YouTubePlayerModule} from '@angular/youtube-player';

// GOF app imports
import { GofHolderComponent } from './gof-holder/gof-holder.component';
import { MapViewComponent } from './gof-holder/map-view/map-view.component';
import { SpeciesdbComponent } from './speciesdb/speciesdb.component';
import { ListspeciesComponent } from './speciesdb/listspecies/listspecies.component';
import { DatabaseService } from './services/database/database.service';
import { FilterspeciePipe } from './pipes/filterspecie.pipe';
import { SpeciesComponent } from './speciesdb/species/species.component';
import { FilterTaxonPipe } from './pipes/filter-taxon.pipe';
import { FilterColloquialPipe } from './pipes/filter-colloquial.pipe';
import { Molecule3DComponent } from './molecule3-d/molecule3-d.component';
import { NcbiinsightsComponent } from './ncbiinsights/ncbiinsights.component';


@NgModule({
  declarations: [
    AppComponent,
    GofHolderComponent,
    HomeComponent,
    FooterComponent,
    NavbarComponent,
    HeaderComponent,
    NotFoundComponent,
    SpeciesdbComponent,
    ListspeciesComponent,
    SpeciesComponent,
    FilterspeciePipe,
    FilterTaxonPipe,
    FilterColloquialPipe,
    Molecule3DComponent,
    NcbiinsightsComponent,
    
  ],
  imports: [
    BrowserModule,
    AppRoutingModule,
    FormsModule,
    HttpClientModule,
    NgxPaginationModule,
    YouTubePlayerModule
  ],
  providers: [DatabaseService],
  bootstrap: [AppComponent]
})
export class AppModule { }
